#!/bin/bash

install_homebrew() {
    echo "Homebrew not found. Installing Homebrew..."
    /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
}

if ! command -v brew &>/dev/null; then
    install_homebrew
else
    echo "Homebrew already installed XD"
fi

echo "Installing cmake git and gsl via homebrew lol"
brew install cmake git gsl

RACIPE_VERSION="1.0"
RACIPE_REPO="https://github.com/simonhb1990/RACIPE-1.0.git"
RACIPE_DIR="RACIPE-$RACIPE_VERSION"

# Clone the RACIPE repository
if [ ! -d "$RACIPE_DIR" ]; then
    echo "Cloning RACIPE repository..."
    git clone $RACIPE_REPO $RACIPE_DIR
else
    echo "RACIPE directory already exists."
fi

cd $RACIPE_DIR

echo "Building RACIPE..."
make

if [ -f "./RACIPE" ]; then
    echo "RACIPE built successfully."
else
    echo "Error: RACIPE build failed."
    exit 1
fi

echo "Adding RACIPE to PATH..."
SHELL_CONFIG=""
if [[ "$SHELL" == */zsh ]]; then
    SHELL_CONFIG="$HOME/.zshrc"
elif [[ "$SHELL" == */bash ]]; then
    SHELL_CONFIG="$HOME/.bash_profile"
fi

if [ -n "$SHELL_CONFIG" ]; then
    if ! grep -q "export PATH=.*$PWD" "$SHELL_CONFIG"; then
        echo "export PATH=\$PATH:$PWD" >>"$SHELL_CONFIG"
        source "$SHELL_CONFIG"
        echo "RACIPE added to PATH in $SHELL_CONFIG."
    else
        echo "RACIPE already in PATH."
    fi
else
    echo "Unknown shell. add RACIPE to PATH manually: export PATH=\$PATH:$PWD"
fi

echo "Verifying RACIPE installation..."
if command -v RACIPE &>/dev/null; then
    echo "RACIPE installed successfully. You can run it using 'RACIPE'."
else
    echo "Error: RACIPE not found in PATH. Please check your PATH settings."
fi
