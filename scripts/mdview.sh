#!/bin/bash
# mdview - View markdown files in browser

if [ -z "$1" ]; then
    echo "Usage: $0 <markdown_file> [browser]"
    echo "Example: $0 README.md"
    echo "         $0 AGENTS.md chromium"
    exit 1
fi

MD_FILE="$1"

# Detect available browser
if command -v google-chrome &> /dev/null; then
    BROWSER="google-chrome"
elif command -v chromium &> /dev/null; then
    BROWSER="chromium"
elif command -v chromium-browser &> /dev/null; then
    BROWSER="chromium-browser"
elif command -v firefox &> /dev/null; then
    BROWSER="firefox"
elif command -v epiphany-browser &> /dev/null; then
    BROWSER="epiphany-browser"
else
    echo "No browser found. Please install one or specify as argument."
    exit 1
fi

# Allow override from argument
if [ -n "$2" ]; then
    BROWSER="$2"
fi

if [ ! -f "$MD_FILE" ]; then
    echo "Error: File '$MD_FILE' not found"
    exit 1
fi

OUTPUT_FILE="$(basename "$MD_FILE" .md).html"

pandoc "$MD_FILE" -o "$OUTPUT_FILE"

echo "Generated: $OUTPUT_FILE"

case "$BROWSER" in
    firefox|chromium|google-chrome|brave|lynx)
        "$BROWSER" "$OUTPUT_FILE" &
        ;;
    *)
        echo "Unknown browser: $BROWSER"
        echo "Using xdg-open instead..."
        xdg-open "$OUTPUT_FILE" &
        ;;
esac

echo "Opening in $BROWSER..."
