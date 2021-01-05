# PREDAIP
PREDAIP uses machine learning algorithms to accurately predict anti-inflammatory peptides of information from peptide sequences.

# Installation
* Install Python 3.7 (https://www.python.org/downloads/) in Linux and Windows.
* Because the program is written in Python 3.7, python 3.7 with the pip tool must be installed first. PREDAIP uses the following dependencies: numpy, pandas, matplotlib, scipy, sklearn version == 0.20.3. You can install these packages first, by the following commands:
```
pip install numpy
pip install pandas
pip install matplotlib
pip install scipy
pip install sklearn==0.20.3
```
* If you meet an error after inputting above commands in Linux, the specific contents are as follows:
</br>Could not install packages due to an EnvironmentError: [Errno 13] Permission denied: '/usr/local/lib/python3.7/dist-packages/sklearn'
Consider using the '--user' option or check the permissions.
</br>Users can change the commands into:
```
pip install numpy --user
pip install pandas --user
pip install matplotlib --user
pip install scipy --user
pip install scikit-learn --user
```

# Running PREDAIP
Download the ‘.py’ files.
Open cmd in Windows or terminal in Linux, then cd to the PREDAIP-master/codes folder which contains predict.py
</br>**For AIP prediction using our model, run:**
</br>`python predict.py -input [custom predicting data in txt format] -threshold [threshold value] -output [predicting results in csv format]` 
</br>**Example:**
</br>`python predict.py -input ../codes/example.txt -threshold 0.5 -output ../codes/results.csv`
</br>-output is optional parameter, while -input and -threshold are required parameters. Prediction results will show in the cmd or terminal, and if you don't want to save results, you need not input -output.

</br>**Example:**
</br>`python predict.py -input ../codes/example.txt -threshold 0.5`

</br>**For details of -input, -threshold and -output, run:**
</br>`python predict.py -h`

# Announcements
* In order to obtain the prediction results, please save the query peptide sequences and peptide names in txt format. Users can refer to the example.txt under the codes folder. Also of note, each peptide name should be added by '>', otherwise the program will occur error. 
* To save the prediction results, the -ouput should be in csv format.
