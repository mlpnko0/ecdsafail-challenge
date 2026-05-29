#!/usr/bin/env bash
# Ensure Rust + a working C linker are installed. Idempotent: no-op if
# `cargo` and a supported C compiler/linker are already on PATH.
set -euo pipefail

SUDO=""
if [[ ${EUID:-$(id -u)} -ne 0 ]] && command -v sudo >/dev/null 2>&1; then
  SUDO="sudo"
fi

# shellcheck disable=SC1091
. "$HOME/.cargo/env" 2>/dev/null || true

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

install_system_deps() {
  if command -v apt-get >/dev/null 2>&1; then
    export DEBIAN_FRONTEND=noninteractive
    ${SUDO} apt-get update
    ${SUDO} apt-get install -y --no-install-recommends gcc libc6-dev ca-certificates
  elif command -v dnf >/dev/null 2>&1; then
    ${SUDO} dnf install -y gcc glibc-devel ca-certificates
  elif command -v yum >/dev/null 2>&1; then
    ${SUDO} yum install -y gcc glibc-devel ca-certificates
  elif command -v apk >/dev/null 2>&1; then
    ${SUDO} apk add --no-cache gcc musl-dev ca-certificates
  elif command -v pacman >/dev/null 2>&1; then
    ${SUDO} pacman -Sy --noconfirm gcc ca-certificates
  elif command -v zypper >/dev/null 2>&1; then
    ${SUDO} zypper --non-interactive install gcc glibc-devel ca-certificates
  elif command -v brew >/dev/null 2>&1; then
    : # macOS: cc comes from Xcode CLT, which `xcode-select --install` handles.
  else
    return 1
  fi
}

# 1. System deps: a C compiler/linker.
#    Cargo needs a linker, and some crates may shell out through the `cc` crate.
if ! find_c_compiler >/dev/null 2>&1; then
  if ! install_system_deps; then
    cat >&2 <<'EOF'
setup.sh: failed to install system dependencies.

This environment needs a C compiler/linker before Rust can build this
repo. If package metadata cannot be downloaded, enable network/DNS for the
sandbox or use an image that already includes gcc (or clang), glibc-devel,
and ca-certificates.
EOF
    exit 1
  fi
fi

# Prefer a compiler that already exists. `benchmark.sh` repeats this detection
# because exported variables from setup.sh do not persist across separate runs.
compiler="$(find_c_compiler || true)"
if [[ -z "${compiler}" ]]; then
  echo "setup.sh: no C compiler found; install gcc or clang" >&2
  exit 1
fi
export CC="${compiler}"

# 2. Rust toolchain.
if ! command -v cargo >/dev/null 2>&1; then
  if ! curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs \
    | sh -s -- -y --profile minimal; then
    cat >&2 <<'EOF'
setup.sh: failed to install Rust with rustup.

If this sandbox has no outbound network/DNS, use an image that already includes
rustup/cargo and this repo's Rust toolchain.
EOF
    exit 1
  fi
fi

# shellcheck disable=SC1091
. "$HOME/.cargo/env" 2>/dev/null || true

if ! command -v cargo >/dev/null 2>&1; then
  echo "setup.sh: cargo is still not available after setup" >&2
  exit 1
fi

channel="$(pinned_rust_channel)"
if [[ -n "${channel}" ]] && command -v rustup >/dev/null 2>&1; then
  toolchain="$(installed_toolchain_for_channel "${channel}" || true)"
  if [[ -z "${toolchain}" ]]; then
    rustup toolchain install "${channel}" --profile minimal
    toolchain="$(installed_toolchain_for_channel "${channel}" || true)"
  fi

  if [[ -z "${toolchain}" ]]; then
    echo "setup.sh: failed to install Rust toolchain '${channel}'" >&2
    exit 1
  fi

  export RUSTUP_TOOLCHAIN="${toolchain}"
fi

# 3. Populate the Cargo cache and prebuild the benchmark binaries. After this,
#    benchmark.sh should not need network access.
cargo fetch --locked
RUSTFLAGS="-C linker=${compiler}" cargo build --release --locked --bin build_circuit --bin eval_circuit
