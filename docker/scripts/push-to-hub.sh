#!/bin/bash
# Push Docker image to Docker Hub for team sharing

# Configuration - update these for your organization
DOCKER_HUB_USER="${DOCKER_HUB_USER:-nanosystemslab}"
REPO_NAME="${REPO_NAME:-bridge-connection}"
IMAGE_NAME="${1:-me-672:amd64}"
TAG="${2:-latest}"

FULL_NAME="$DOCKER_HUB_USER/$REPO_NAME:$TAG"

echo "========================================="
echo "Docker Hub Image Push"
echo "========================================="
echo "Source image: $IMAGE_NAME"
echo "Target: $FULL_NAME"
echo ""

# Check if logged in to Docker Hub
if ! docker info 2>/dev/null | grep -q "Username"; then
    echo "Please log in to Docker Hub:"
    docker login
fi

# Check if source image exists
if ! docker image inspect "$IMAGE_NAME" >/dev/null 2>&1; then
    echo "Error: Image '$IMAGE_NAME' not found!"
    echo "Available images:"
    docker images
    exit 1
fi

# Tag the image
echo "Tagging image..."
docker tag "$IMAGE_NAME" "$FULL_NAME"

# Push to Docker Hub
echo "Pushing to Docker Hub (this may take several minutes)..."
docker push "$FULL_NAME"

if [ $? -eq 0 ]; then
    echo ""
    echo "✓ Image pushed successfully!"
    echo ""
    echo "Team members can now pull the image with:"
    echo "  docker pull $FULL_NAME"
    echo ""
    echo "Update docker-compose.yml to use:"
    echo "  image: $FULL_NAME"
else
    echo "Error: Failed to push image!"
    echo "Make sure you have access to $DOCKER_HUB_USER organization on Docker Hub"
    exit 1
fi