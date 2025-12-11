#!/bin/bash
# Load Docker image from tar.gz file

INPUT_FILE="${1:-docker/images/bridge-connection.tar.gz}"

echo "Loading Docker image from: $INPUT_FILE"

# Check if file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: File '$INPUT_FILE' not found!"
    echo ""
    echo "Please download the image file first:"
    echo "  1. From team's cloud storage (Google Drive/Dropbox)"
    echo "  2. Or pull from Docker Hub:"
    echo "     docker pull nanosystemslab/bridge-connection:latest"
    exit 1
fi

# Show file info
SIZE=$(du -h "$INPUT_FILE" | cut -f1)
echo "File size: $SIZE"

# Load image
echo "Loading image (this may take a few minutes)..."
docker load < "$INPUT_FILE"

if [ $? -eq 0 ]; then
    echo "✓ Image loaded successfully!"
    echo ""
    echo "Available images:"
    docker images
    echo ""
    echo "To run the container:"
    echo "  docker compose up"
    echo "Or:"
    echo "  docker run --rm -it -p 8888:8888 \\"
    echo "    -v \"\$PWD\":/home/fenics/shared bridge-connection:dolfinx"
else
    echo "Error: Failed to load image!"
    exit 1
fi