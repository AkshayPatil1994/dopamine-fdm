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

# Wiki pages are addressed without the .md extension (a link to Page.md opens the raw file), and links that
# leave docs/ must point at the repository. docs/ itself keeps the .md links so it also works when browsed in the repo.
repo_web="${remote_url%.git}"
repo_web="https://github.com/${repo_web#*github.com[:/]}/blob/main/"
find . -name '*.md' -not -path './.git/*' -print0 | xargs -0 perl -pi -e \
  "s{\]\(([A-Za-z0-9_-]+)\.md(#[^)]*)?\)}{](\$1\$2)}g; s{\]\(\.\./}{](${repo_web}}g"

git add -A
if git diff --cached --quiet; then
  echo "No changes to sync."
  exit 0
fi

git commit -m "Sync wiki from docs/"
git push
echo "Wiki updated."
