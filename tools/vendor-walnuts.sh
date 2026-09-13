#!/bin/sh
# Re-vendor the WALNUTS headers from upstream.
#
# Usage:
#   tools/vendor-walnuts.sh <upstream-sha> [destdir]
#
# Fetches include/walnutpie/ from upstream at <upstream-sha>, trims it to the
# file list in inst/include/WALNUTS_VERSION, rewrites the include paths from
# walnutpie/ to walnuts/, regenerates the umbrella walnuts.hpp, refreshes
# WALNUTS_LICENSE, reapplies every patch in tools/walnuts-patches/, and writes
# an updated WALNUTS_VERSION with the new sha and dates.
#
# With no destdir it writes inst/include in place. With one it writes a
# self-contained copy (walnuts/, walnuts.hpp, WALNUTS_LICENSE, WALNUTS_VERSION)
# there instead, which is how the result is diffed against what is checked in:
#
#   tools/vendor-walnuts.sh f3c1833 /tmp/v && diff -ru inst/include /tmp/v
#
# The patches are the only hand edits carried across updates. Everything else
# under inst/include/walnuts/ is upstream verbatim apart from the include-path
# rewrite, and a diff that shows anything else means this script is stale.
#
# Needs curl and tar. Uses the GitHub API for the commit date; set
# WALNUTS_COMMIT_DATE=YYYY-MM-DD to supply it when the API is unreachable.

set -eu

REPO=flatironinstitute/walnuts

usage() {
  echo "usage: tools/vendor-walnuts.sh <upstream-sha> [destdir]" >&2
  exit 2
}

[ $# -ge 1 ] && [ $# -le 2 ] || usage
sha=$1

root=$(cd "$(dirname "$0")/.." && pwd)
version_file="$root/inst/include/WALNUTS_VERSION"
patch_dir="$root/tools/walnuts-patches"
dest=${2:-"$root/inst/include"}

[ -f "$version_file" ] || { echo "missing $version_file" >&2; exit 1; }

# The vendored set, one file per line between the markers in WALNUTS_VERSION.
files=$(awk '/BEGIN VENDORED FILES/ { on = 1; next }
             /END VENDORED FILES/   { on = 0 }
             on                     { gsub(/[ \t]/, ""); if ($0 != "") print }' \
        "$version_file")
[ -n "$files" ] || { echo "no vendored file list in $version_file" >&2; exit 1; }

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT INT TERM

echo "fetching $REPO at $sha"
curl -sfL "https://codeload.github.com/$REPO/tar.gz/$sha" | tar -xzf - -C "$tmp"
up=$(echo "$tmp"/*/include/walnutpie)
[ -d "$up" ] || { echo "no include/walnutpie in the fetched tree" >&2; exit 1; }

# Provenance. The API reports both dates in UTC; the second is the committer's.
if [ -n "${WALNUTS_COMMIT_DATE:-}" ]; then
  full_sha=$sha
  commit_date=$WALNUTS_COMMIT_DATE
else
  meta=$(curl -sfL "https://api.github.com/repos/$REPO/commits/$sha")
  full_sha=$(printf '%s\n' "$meta" |
             sed -n 's/.*"sha": "\([0-9a-f]\{40\}\)".*/\1/p' | sed -n 1p)
  commit_date=$(printf '%s\n' "$meta" |
                sed -n 's/.*"date": "\([0-9][0-9-]*\)T.*/\1/p' | sed -n 2p)
  if [ -z "$full_sha" ] || [ -z "$commit_date" ]; then
    echo "could not read the sha and date from the GitHub API;" >&2
    echo "rerun with WALNUTS_COMMIT_DATE=YYYY-MM-DD and a full sha" >&2
    exit 1
  fi
fi

mkdir -p "$dest/walnuts"

# Trim and rewrite the include paths.
for f in $files; do
  [ -f "$up/$f" ] || { echo "upstream has no $f at $sha" >&2; exit 1; }
  sed 's|include "walnutpie/|include "walnuts/|' "$up/$f" > "$dest/walnuts/$f"
done

# The umbrella, upstream's own with the dropped headers filtered out.
awk -v keep="$(echo $files)" '
  BEGIN { n = split(keep, a, " "); for (i = 1; i <= n; i++) want[a[i]] = 1 }
  /^#include "walnutpie\// {
    h = $0; sub(/^#include "walnutpie\//, "", h); sub(/".*$/, "", h)
    if (!(h in want)) next
    sub(/walnutpie\//, "walnuts/")
  }
  { print }' "$up.hpp" > "$dest/walnuts.hpp"

cp "$(dirname "$up")/../LICENSE" "$dest/WALNUTS_LICENSE"

# Reapply the local patches, in name order.
for p in "$patch_dir"/*.patch; do
  [ -f "$p" ] || continue
  echo "applying $(basename "$p")"
  if command -v git > /dev/null 2>&1; then
    (cd "$dest" && git apply -p1 "$p")
  else
    (cd "$dest" && patch -p1 --quiet < "$p")
  fi
done

# Provenance lines. Written to the destination, so a scratch run produces a
# WALNUTS_VERSION to diff rather than editing the tree's.
sed -e "s|^Commit:    .*|Commit:    $full_sha|" \
    -e "s|^Date:      .*|Date:      $commit_date (upstream commit date, UTC)|" \
    -e "s|^Vendored:  .*|Vendored:  $(date -u +%Y-%m-%d)|" \
    "$version_file" > "$dest/WALNUTS_VERSION.new"
mv "$dest/WALNUTS_VERSION.new" "$dest/WALNUTS_VERSION"

echo "vendored $REPO@$full_sha ($commit_date) into $dest"
echo "the notes in WALNUTS_VERSION describe the PREVIOUS upstream; reread them"
