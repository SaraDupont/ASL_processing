
sudo python get-pip.py

sudo pip install --requirement requirements.txt --ignore-installed

if ! type fsl > /dev/null; then
	python fslinstaller.py
fi

/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

brew install dcm2niix

echo 'Done ! Close this window and open a new terminal.'