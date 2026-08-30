#!/usr/bin/env bash
# Sync docs/ to the GitHub wiki. Requires SSH push access to the repo.
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
remote_url="$(git -C "$repo_root" remote get-url origin)"
wiki_url="${remote_url%.git}.wiki.git"

work_dir="$(mktemp -d)"
trap 'rm -rf "$work_dir"' EXIT

git clone "$wiki_url" "$work_dir/wiki"
rsync -av --delete --exclude '.git' "$repo_root/docs/" "$work_dir/wiki/"

cd "$work_dir/wiki"
git add -A
if git diff --cached --quiet; then
  echo "No changes to sync."
  exit 0
fi

git commit -m "Sync wiki from docs/"
git push
echo "Wiki updated."
