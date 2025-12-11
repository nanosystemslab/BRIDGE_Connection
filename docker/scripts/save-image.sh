#!/bin/bash
# Save Docker image to tar.gz file for sharing

IMAGE_NAME="${1:-me-672:amd64}"
OUTPUT_FILE="${2:-docker/images/me-672_amd64.tar.gz}"

echo "Saving Docker image: $IMAGE_NAME"
echo "Output file: $OUTPUT_FILE"

# Check if image exists
if ! docker image inspect "$IMAGE_NAME" >/dev/null 2>&1; then
    echo "Error: Image '$IMAGE_NAME' not found!"
    echo "Available images:"
    docker images
    exit 1
fi

# Create directory if it doesn't exist
mkdir -p "$(dirname "$OUTPUT_FILE")"

# Save and compress image
echo "Saving image (this may take a few minutes)..."
docker save "$IMAGE_NAME" | gzip > "$OUTPUT_FILE"

# Check file size
if [ -f "$OUTPUT_FILE" ]; then
    SIZE=$(du -h "$OUTPUT_FILE" | cut -f1)
    echo "✓ Image saved successfully!"
    echo "  File: $OUTPUT_FILE"
    echo "  Size: $SIZE"
    echo ""
    echo "⚠️  WARNING: This file is too large for Git!"
    echo "   Upload to cloud storage or Docker Hub instead."
    echo ""
    echo "To share via Docker Hub:"
    echo "  docker tag $IMAGE_NAME yourusername/bridge-connection:latest"
    echo "  docker push yourusername/bridge-connection:latest"
else
    echo "Error: Failed to save image!"
    exit 1
fi