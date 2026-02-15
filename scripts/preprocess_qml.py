import os
import sys
import tarfile
import urllib.request
import numpy as np
from PIL import Image
import shutil
import io

# Config
DATASET_BASE_URL = "https://github.com/sueiras/handwritting_characters_database/raw/master/"
PARTS = ["curated.tar.gz.01", "curated.tar.gz.02"]
TEMP_DIR = "temp_dataset"
CENTROIDS_FILE = "centroids.bin"
TEST_SAMPLES_FILE = "test_samples.bin"
IMG_WIDTH = 4
IMG_HEIGHT = 5
NUM_FEATURES = IMG_WIDTH * IMG_HEIGHT
START_ASCII = 33
END_ASCII = 126

def download_dataset():
    if not os.path.exists("curated.tar.gz"):
        print(f"Downloading dataset parts from {DATASET_BASE_URL}...")
        try:
            # Download parts
            for part in PARTS:
                if not os.path.exists(part):
                    url = DATASET_BASE_URL + part
                    print(f"Downloading {part}...")
                    urllib.request.urlretrieve(url, part)

            # Concatenate
            print("Concatenating parts...")
            with open("curated.tar.gz", "wb") as outfile:
                for part in PARTS:
                    with open(part, "rb") as infile:
                        outfile.write(infile.read())
            print("Download and concatenation complete.")

        except Exception as e:
            print(f"Error downloading: {e}")
            sys.exit(1)
    else:
        print("Dataset already downloaded.")

def extract_dataset():
    if not os.path.exists(TEMP_DIR):
        print("Extracting dataset...")
        os.makedirs(TEMP_DIR, exist_ok=True)
        try:
            with tarfile.open("curated.tar.gz", "r:gz") as tar:
                tar.extractall(path=TEMP_DIR)
            print("Extraction complete.")
        except Exception as e:
            print(f"Error extracting: {e}")
            sys.exit(1)

def process_image(img_path):
    try:
        with Image.open(img_path) as img:
            # 1. Convert to grayscale
            img = img.convert('L')

            # 2. Invert (assuming dark text on light background -> light text on dark background)
            # The dataset examples (MNIST style) are usually white on black.
            # If these are real handwriting on paper, they are dark on light.
            # Let's inspect a pixel or assume standard inversion is needed for "feature=1"
            # PIL ImageOps.invert requires RGB or L.
            from PIL import ImageOps
            img = ImageOps.invert(img)

            # 3. Thresholding (Binarize)
            # Pixels < 128 become 0, >= 128 become 255
            threshold = 128
            img = img.point(lambda p: 255 if p > threshold else 0)

            # 4. Center of Mass (Bounding Box + Centering)
            # Get bounding box of non-zero regions
            bbox = img.getbbox()
            if bbox:
                # Crop to content
                img = img.crop(bbox)

                # Create a new square image to place content in center
                # Make it slightly larger to preserve aspect ratio
                w, h = img.size
                max_dim = max(w, h)
                new_size = int(max_dim * 1.2) # 20% padding
                new_img = Image.new('L', (new_size, new_size), 0)

                # Paste centered
                offset_x = (new_size - w) // 2
                offset_y = (new_size - h) // 2
                new_img.paste(img, (offset_x, offset_y))
                img = new_img

            # 5. Resize to Target (4x5)
            # Use High Quality downsampling
            img = img.resize((IMG_WIDTH, IMG_HEIGHT), Image.Resampling.LANCZOS)

            # 6. Normalize to 0-1
            arr = np.array(img, dtype=np.float32) / 255.0

            # 7. Final Threshold to ensure clean 0 or 1 specifics (optional but good for encoding)
            # arr[arr < 0.1] = 0
            # arr[arr > 0.9] = 1

            return arr.flatten()
    except Exception as e:
        print(f"Error processing {img_path}: {e}")
        return None

def main():
    download_dataset()
    extract_dataset()

    centroids = []
    test_samples = []

    # Path to the data inside extracted folder
    # Usually curated.tar.gz extracts to a folder named 'curated' or similar.
    # Let's find it.
    base_path = os.path.join(TEMP_DIR, "curated")
    if not os.path.exists(base_path):
        # Fallback: check if it extracted directly to root or another name
        # Listing dirs
        subdirs = [d for d in os.listdir(TEMP_DIR) if os.path.isdir(os.path.join(TEMP_DIR, d))]
        if len(subdirs) == 1:
            base_path = os.path.join(TEMP_DIR, subdirs[0])
        else:
            print(f"Structure unknown. Dirs: {subdirs}")
            # Try to find number folders directly
            base_path = TEMP_DIR

    print(f"Processing classes {START_ASCII} to {END_ASCII}...")

    count_classes = 0
    final_centroids = {} # Map char id -> vector
    all_train_samples = []

    for ascii_code in range(START_ASCII, END_ASCII + 1):
        folder_name = str(ascii_code)
        folder_path = os.path.join(base_path, folder_name)

        if not os.path.isdir(folder_path):
            # print(f"Warning: Class {folder_name} ({chr(ascii_code)}) not found.")
            # We must supply a zero vector or skip?
            # Better to store a dummy to keep indexing aligned if we used arrays,
            # but we can store ID in the binary file.
            continue

        images = []
        for f in os.listdir(folder_path):
            if f.startswith("._"):
                try:
                    os.remove(os.path.join(folder_path, f))
                except OSError:
                    pass
                continue
            if f.endswith('.png'):
                images.append(os.path.join(folder_path, f))

        if not images:
            continue

        # Split Train/Test (Simple 90/10 split)
        # Or just use all for centroid and save a few specific ones for testing
        vectors = []
        for img_path in images:
            vec = process_image(img_path)
            if vec is not None:
                vectors.append(vec)

        if not vectors:
            continue

        vectors = np.array(vectors)

        # Calculate Centroid
        centroid = np.mean(vectors, axis=0)
        final_centroids[ascii_code] = centroid

        # Collect k-NN Training Data (Exclude test samples)
        num_test = min(5, len(vectors))
        train_vectors = vectors[:-num_test] if num_test > 0 else vectors

        # Limit to 50 training samples per class
        if len(train_vectors) > 50:
            train_vectors = train_vectors[:50]

        for v in train_vectors:
            all_train_samples.append((ascii_code, v))

        # Save last 5 as test samples
        num_test = min(5, len(vectors))
        for i in range(num_test):
            # Format: ID (1 byte) + Vector (20 floats)
            test_samples.append((ascii_code, vectors[-(i+1)]))

        count_classes += 1
        print(f"Processed '{chr(ascii_code)}' ({ascii_code}): {len(vectors)} images.")

    print(f"\nFinal Statistics: {count_classes} classes processed.")

    # Write Full Training Set
    # Format: NumSamples (int), [ClassID (int), Vector (20 floats)]...
    TRAIN_FILE = "train_samples.bin"
    with open(TRAIN_FILE, "wb") as f:
        # Write total count
        f.write(np.array([len(all_train_samples)], dtype=np.int32).tobytes())
        for class_id, vec in all_train_samples:
            f.write(np.array([class_id], dtype=np.int32).tobytes())
            f.write(vec.astype(np.float32).tobytes())
    print(f"Saved {TRAIN_FILE} ({len(all_train_samples)} samples)")

    # Write Test Samples File
    # Format: NumSamples (int), [ClassID (int), Vector (20 floats)]...
    with open(TEST_SAMPLES_FILE, "wb") as f:
        f.write(np.array([len(test_samples)], dtype=np.int32).tobytes())
        for ascii_code, vec in test_samples:
            f.write(np.array([ascii_code], dtype=np.int32).tobytes())
            f.write(vec.astype(np.float32).tobytes())
    print(f"Saved {TEST_SAMPLES_FILE} ({len(test_samples)} samples)")

    # Cleanup intermediate files
    print("Cleaning up intermediate files...")
    if os.path.exists(TEMP_DIR):
        shutil.rmtree(TEMP_DIR)

    # Remove tar files
    for part in PARTS:
        if os.path.exists(part):
            os.remove(part)
    if os.path.exists("curated.tar.gz"):
        os.remove("curated.tar.gz")
    print("Cleanup complete.")

if __name__ == "__main__":
    main()
