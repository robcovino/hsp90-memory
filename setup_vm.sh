#!/usr/bin/env bash
# One-shot setup of a fresh Ubuntu VM for the hsp90-memory project.
# Run as the normal user (not root):   bash setup_vm.sh
# Idempotent: safe to run again. Takes a few minutes the first time.
set -euo pipefail

REPO_URL="${REPO_URL:-https://github.com/robcovino/hsp90-memory.git}"
SCRIPT_URL="https://raw.githubusercontent.com/robcovino/hsp90-memory/main/setup_vm.sh"

# Commercial clouds (DigitalOcean, Hetzner, ...) log you in as root. Create a normal user 'ubuntu' with passwordless
# sudo, copy the SSH key over, and continue as that user, so the rest is identical to de.NBI (which uses 'ubuntu').
if [ "$(id -u)" = 0 ]; then
  echo "== running as root: creating user 'ubuntu' and continuing as that user"
  id ubuntu > /dev/null 2>&1 || adduser --disabled-password --gecos "" ubuntu
  usermod -aG sudo ubuntu
  echo 'ubuntu ALL=(ALL) NOPASSWD:ALL' > /etc/sudoers.d/90-ubuntu && chmod 440 /etc/sudoers.d/90-ubuntu
  mkdir -p /home/ubuntu/.ssh
  [ -f /root/.ssh/authorized_keys ] && cp -n /root/.ssh/authorized_keys /home/ubuntu/.ssh/authorized_keys
  chown -R ubuntu:ubuntu /home/ubuntu/.ssh && chmod 700 /home/ubuntu/.ssh
  echo "   (from now on log in as:  ssh ubuntu@<this VM's address>)"
  exec su - ubuntu -c "curl -fsSL $SCRIPT_URL | bash"
fi

DIR="$HOME/hsp90-memory"

echo "== system packages"
sudo apt-get update -qq
sudo apt-get install -y -qq git tmux curl wget python3 python3-venv python3-pip build-essential > /dev/null

echo "== Claude Code (native installer)"
if ! command -v claude > /dev/null 2>&1 && [ ! -x "$HOME/.local/bin/claude" ]; then
  curl -fsSL https://claude.ai/install.sh | bash
fi
export PATH="$HOME/.local/bin:$PATH"
grep -q '.local/bin' "$HOME/.bashrc" || echo 'export PATH="$HOME/.local/bin:$PATH"' >> "$HOME/.bashrc"

echo "== repository"
if [ -d "$DIR/.git" ]; then
  git -C "$DIR" pull --ff-only
else
  git clone "$REPO_URL" "$DIR"
fi
cd "$DIR"

echo "== python environment"
[ -d .venv ] || python3 -m venv .venv
# shellcheck disable=SC1091
source .venv/bin/activate
pip install -q --upgrade pip
pip install -q -e ".[dev,sbi]" numba

echo "== data (Zenodo 12697153)"
mkdir -p data
for i in 1 2 3 4; do
  f="data/${i}Thermal.txt"
  [ -s "$f" ] || wget -q "https://zenodo.org/records/12697153/files/${i}Thermal.txt?download=1" -O "$f"
done
ls -l data/*.txt

echo "== tests"
python -m pytest -q

cat <<'EOF'

Setup complete. Next, two commands:

  1) sign in once (as user ubuntu): cd ~/hsp90-memory && source .venv/bin/activate && claude
     (type /login, open the printed URL on your laptop, paste the code back; answer "yes" to trusting the folder; then type /exit)

  2) hand the machine to Remote Control:   bash ~/hsp90-memory/start_remote.sh

Then open https://claude.ai/code on any device: the session "hsp90-memory" is there.
EOF
