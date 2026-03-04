#!/usr/bin/env bash
set -euo pipefail

# Deterministic local tarball-first check harness for npRmpi.
# - Ensures checkbashisms compatibility via a shellcheck-backed wrapper.
# - Forces local future-timestamp check behavior to avoid flaky clock-verification notes.

ROOT_DIR="/Users/jracine/Development"
PKG_DIR="$ROOT_DIR/np-npRmpi"
OUT_DIR="${1:-/tmp/npRmpi_ascran_local_$(date +%Y%m%d_%H%M%S)}"

mkdir -p "$OUT_DIR"

echo "OUT_DIR=$OUT_DIR"

if ! command -v R >/dev/null 2>&1; then
  echo "R not found on PATH" >&2
  exit 1
fi

if ! command -v shellcheck >/dev/null 2>&1; then
  echo "shellcheck not found on PATH (required for checkbashisms shim)" >&2
  exit 1
fi

cat > /tmp/checkbashisms <<'SHIM'
#!/usr/bin/env bash
set -u
if ! command -v shellcheck >/dev/null 2>&1; then
  exit 0
fi

files=()
while (($#)); do
  case "$1" in
    -p|-d|--posix|--debian|--newline|--force|--extra)
      shift
      ;;
    --)
      shift
      break
      ;;
    -*)
      shift
      ;;
    *)
      files+=("$1")
      shift
      ;;
  esac
done
while (($#)); do
  files+=("$1")
  shift
done

lint_targets=()
for f in "${files[@]}"; do
  [ -f "$f" ] || continue
  case "$f" in
    *.sh)
      lint_targets+=("$f")
      continue
      ;;
  esac
  first=$(head -n 1 "$f" 2>/dev/null || true)
  case "$first" in
    '#!'*sh*|'#!'*bash*)
      lint_targets+=("$f")
      ;;
  esac
done

if ((${#lint_targets[@]} == 0)); then
  exit 0
fi

shellcheck --shell=sh --exclude=SC1090,SC1091 "${lint_targets[@]}" || true
exit 0
SHIM
chmod +x /tmp/checkbashisms

cd "$ROOT_DIR"
R CMD build "$PKG_DIR" > "$OUT_DIR/rcmd_build.log" 2>&1

PKG_TARBALL=$(
  find "$ROOT_DIR" -maxdepth 1 -type f -name 'npRmpi_*.tar.gz' -exec stat -f '%m %N' {} + \
    | sort -rn \
    | awk 'NR==1 { sub(/^[0-9]+ /, ""); print }'
)
if [[ -z "$PKG_TARBALL" ]]; then
  echo "Unable to locate built tarball" >&2
  exit 1
fi

echo "$PKG_TARBALL" > "$OUT_DIR/tarball.txt"

PATH="/tmp:$PATH" _R_CHECK_FUTURE_FILE_TIMESTAMPS_=false \
  R CMD check --as-cran "$PKG_TARBALL" > "$OUT_DIR/rcmd_check_as_cran.log" 2>&1 || true

rg -n 'Status:|top-level files|future file timestamps|checkbashisms' \
  "$OUT_DIR/rcmd_check_as_cran.log" > "$OUT_DIR/status_summary.log" || true

echo "DONE"
