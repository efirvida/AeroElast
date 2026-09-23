#!/bin/bash
# Run the WINE_AppImage (mmtrt) wine binaries on hosts with old glibc (RHEL 8).
# Bypasses the AppImage AppRun runtime, which breaks on RHEL 8's backported
# glibc (_dl_call_libc_early_init assertion), by invoking the binaries
# directly with the AppImage's bundled compat glibc 2.31 loader.
set -u
APP_DIR=/scratch/leahk/eduardo.donestevez/AppDir
export WINEPREFIX="${WINEPREFIX:-$HOME/.wine-appimage-stable}"
export WINEDLLPATH="$APP_DIR/opt/wine-stable/lib/wine/x86_64-unix:$APP_DIR/opt/wine-stable/lib/wine"
export WINEDEBUG="${WINEDEBUG:--all}"

if [ "$#" -lt 1 ]; then
  echo "usage: $0 <wine-tool> [args...]" >&2
  echo "e.g.:   $0 wine /path/to/program.exe" >&2
  exit 1
fi

BIN="$APP_DIR/opt/wine-stable/bin/$1"
shift
exec "$APP_DIR/runtime/compat/lib64/ld-linux-x86-64.so.2" \
  --library-path "$APP_DIR/runtime/compat/lib64:$APP_DIR/runtime/compat/lib/x86_64-linux-gnu:$APP_DIR/lib/x86_64-linux-gnu:$APP_DIR/usr/lib/x86_64-linux-gnu:$APP_DIR/opt/wine-stable/lib/wine/x86_64-unix:$APP_DIR/opt/wine-stable/lib/wine" \
  "$BIN" "$@"