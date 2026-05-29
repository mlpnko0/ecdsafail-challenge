#!/usr/bin/env bash
# Benchmark a submission.
#
#   1. Wipe stale ops.bin / score.json so a contestant cannot pre-seed them.
#   2. Run build_circuit (UNTRUSTED — runs contestant code in src/point_add)
#      in its own process group; on exit, kill the whole group so a forked
#      child cannot survive and overwrite score.json after step 4.
#   3. Verify ops.bin exists. If build_circuit exited early or crashed,
#      the file is missing and we fail closed.
#   4. Run eval_circuit (TRUSTED — never imports contestant code) which
#      re-simulates the op stream, validates correctness/reversibility/
#      phase, and writes the canonical score.json.
#
# All command-line arguments are forwarded to eval_circuit (e.g. --note ...).
set -euo pipefail

# shellcheck disable=SC1091
. "$HOME/.cargo/env" 2>/dev/null || true

pinned_rust_channel() {
  local channel=""
  if [[ -f rust-toolchain ]]; then
    channel="$(sed -n 's/^[[:space:]]*channel[[:space:]]*=[[:space:]]*"\([^"]*\)".*/\1/p' rust-toolchain | sed -n '1p')"
    if [[ -z "${channel}" ]]; then
      channel="$(sed -n '1s/^[[:space:]]*\([^[:space:]]*\)[[:space:]]*$/\1/p' rust-toolchain)"
    fi
  fi
  printf '%s\n' "${channel}"
}

installed_toolchain_for_channel() {
  local channel="$1"
  local line toolchain

  [[ -n "${channel}" ]] || return 1
  command -v rustup >/dev/null 2>&1 || return 1

  while IFS= read -r line; do
    toolchain="${line%% *}"
    if [[ "${toolchain}" == "${channel}" || "${toolchain}" == "${channel}-"* ]]; then
      printf '%s\n' "${toolchain}"
      return 0
    fi
  done < <(rustup toolchain list 2>/dev/null)

  return 1
}

require_offline_rust_toolchain() {
  if [[ -n "${RUSTUP_TOOLCHAIN:-}" ]] || ! command -v rustup >/dev/null 2>&1; then
    return 0
  fi

  local channel toolchain
  channel="$(pinned_rust_channel)"
  [[ -n "${channel}" ]] || return 0

  toolchain="$(installed_toolchain_for_channel "${channel}" || true)"
  if [[ -z "${toolchain}" ]]; then
    echo "!! pinned Rust toolchain '${channel}' is not installed; run ./setup.sh before offline benchmarking" >&2
    exit 1
  fi

  export RUSTUP_TOOLCHAIN="${toolchain}"
}

find_c_compiler() {
  if [[ -n "${CC:-}" ]] && command -v "${CC}" >/dev/null 2>&1; then
    command -v "${CC}"
    return 0
  fi

  local candidate
  for candidate in gcc cc clang; do
    if command -v "${candidate}" >/dev/null 2>&1; then
      command -v "${candidate}"
      return 0
    fi
  done

  return 1
}

compiler="$(find_c_compiler || true)"
if [[ -z "${compiler}" ]]; then
  echo "!! no C compiler/linker found; run ./setup.sh or install gcc/clang" >&2
  exit 1
fi
export CC="${compiler}"
require_offline_rust_toolchain

if ! command -v cargo >/dev/null 2>&1; then
  echo "!! cargo not found; run ./setup.sh before offline benchmarking" >&2
  exit 1
fi

export CARGO_NET_OFFLINE=true

# 1. Clean slate.
rm -f ops.bin score.json

# Make sure both binaries are present (cheap rebuild — no-op if up to date).
RUSTFLAGS="-C linker=${compiler}" cargo build --release --locked --offline --bin build_circuit --bin eval_circuit

# 2. Run build_circuit in its own process group, then nuke the group.
#    `setsid` puts it in a fresh pgid; `kill -KILL -<pgid>` reaches every
#    child including double-forked daemons. We trap so we always reap.
#    On macOS `setsid` is missing — fall back to `set -m` + `kill %1`.
cleanup_pgid=""
cleanup() {
  if [[ -n "${cleanup_pgid}" ]]; then
    kill -KILL -"${cleanup_pgid}" 2>/dev/null || true
  fi
}
trap cleanup EXIT

if command -v setsid >/dev/null 2>&1; then
  setsid ./target/release/build_circuit &
  build_pid=$!
  cleanup_pgid="${build_pid}"
  set +e
  wait "${build_pid}"
  build_status=$?
  set -e
  kill -KILL -"${cleanup_pgid}" 2>/dev/null || true
  cleanup_pgid=""
else
  # Fallback: bash job control puts the background pipeline in its own pgid.
  set -m
  ./target/release/build_circuit &
  build_pid=$!
  cleanup_pgid="${build_pid}"
  set +e
  wait "${build_pid}"
  build_status=$?
  set -e
  kill -KILL -"${cleanup_pgid}" 2>/dev/null || true
  cleanup_pgid=""
  set +m
fi

if [[ "${build_status}" -ne 0 ]]; then
  echo "!! build_circuit exited with status ${build_status}" >&2
  exit "${build_status}"
fi

# 3. Verify ops.bin actually got produced.
if [[ ! -s ops.bin ]]; then
  echo "!! build_circuit did not produce ops.bin" >&2
  exit 1
fi

# 4. Trusted scoring stage.
./target/release/eval_circuit "$@"
