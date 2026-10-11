#!/usr/bin/env bash
# Start (or re-attach to) the Claude Code Remote Control server for this repository inside tmux,
# so it keeps running after you log out. Run after setup_vm.sh and a one-time `claude` /login.
set -euo pipefail
export PATH="$HOME/.local/bin:$PATH"
DIR="$HOME/hsp90-memory"
if tmux has-session -t cc 2> /dev/null; then
  echo "Remote Control is already running in tmux session 'cc'. Attach with: tmux attach -t cc   (detach: Ctrl-b then d)"
  exit 0
fi
tmux new-session -d -s cc -c "$DIR" \
  "source .venv/bin/activate && claude remote-control --name hsp90-memory --permission-mode acceptEdits; exec bash"
echo "Started. The session appears at https://claude.ai/code within a minute (and in the Claude phone app)."
echo "To watch it in the terminal: tmux attach -t cc   (detach without stopping: Ctrl-b then d)"
