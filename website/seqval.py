import sys, os, shutil, subprocess, shlex, time
from smb.SMBConnection import SMBConnection
from smb.base import NotConnectedError, NotReadyError, SMBTimeout, SharedFile
from flask import Flask, render_template, request, jsonify
import tempfile

# Initialize the Flask application
app = Flask(__name__)

#SMB credentials - SECRET
username = "secret"
password = "secret"
#Specific server information
myRequestIdentifier = "miseqvalpipeline"
serverName = "SMB"
domain = ""
host = "secret.jbei.org"
port = 139
sharedFolder = "miseq"

class idtext():
	id = 0
	text = ""

	def __init__(self, id, text): 
		self.id = id
		self.text = text
	
	def serialize(self):
		return {
		'id': self.id, 
		'text': self.text,
		}

#
#Return Main Libraries
#

@app.route('/return_mainlibraries', methods=['GET', 'POST'])
def return_mainlibraries():
	mainLibraryNames = []
	query = request.form['query']
	try:
		conn = SMBConnection(username, password, myRequestIdentifier, serverName, domain=domain, use_ntlm_v2 = True)
		conn.connect(host, port)

		mainLibraryFolders = conn.listPath(sharedFolder, '/MiSeqOutput')
		for f in mainLibraryFolders:
			#Ignore .DS_Store
			folderName = f.filename
			if folderName.startswith("."):
				continue
			if query.lower() in folderName.lower():
				mainLibraryNames.append(idtext(folderName, folderName))
	except Exception as ex:
		return jsonify(result=str(ex))
	return jsonify(result=[e.serialize() for e in mainLibraryNames])


#
#Return Sub Libraries
#

@app.route('/return_sampleIDs', methods=['GET', 'POST'])
def return_sampleIDs():
	sampleIDs = []
	query = request.form['query']
	mainLibraryFolder = request.form['mainLibraryFolder']
	try:
		conn = SMBConnection(username, password, myRequestIdentifier, serverName, domain=domain, use_ntlm_v2 = True)
		conn.connect(host, port)

		sampleSheetCSV = tempfile.NamedTemporaryFile()
		pathTo = 'MiSeqOutput/'+mainLibraryFolder+'/SampleSheet.csv'
		sampleSheetCSV_attributes, sampleSheetCSV_size = conn.retrieveFile(sharedFolder, pathTo, sampleSheetCSV)

		sampleSheetCSV.seek(0)

		fileContents = sampleSheetCSV.read()
		uniqueLines = fileContents.replace("\r\n", '\n').replace("\r", '\n').split("\n")

		counter = 0
		for line in uniqueLines:
			#sampleIDs.append(idtext(line, line))
			if (line.startswith("[Data]") or counter==1):
				counter+=1
				continue
			#Two lines after [Data] line, first sampleIDs is encountered
			if (counter==2):
				sampleID = line[:line.find(",")]
				if (query.lower() in sampleID.lower()) and not sampleID=="": #Not blank line
					sampleIDs.append(idtext(sampleID, sampleID))
	except Exception as ex:
		exc_type, exc_obj, exc_tb = sys.exc_info()
		fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
		return jsonify(result=(exc_type, fname, exc_tb.tb_lineno))
	return jsonify(result=[e.serialize() for e in sampleIDs])


#
#Return ICE reference entries
#


@app.route('/return_iceEntries', methods=['GET', 'POST'])
def return_iceEntries():
	iceEntries = []
	return jsonify(result=[e.serialize() for e in iceEntries])

#
#Run Pipeline.py
#


@app.route('/run_pipeline', methods=['GET', 'POST'])
def run_pipeline():
	email = request.form['email']
	mainLibraryFolder = request.form['mainLibrary']
	subLibraryIDs = request.form['sampleIDs'].replace(",", " ")
	referenceSequences = 'references' #Uncomment this to get POSTed form data: request.form['referenceSequence']

	error = False
	emailError = False
	mainLibraryFolderError = False
	subLibraryIDsError = False
	referenceSequencesError = False

	if (email==""):
		error = True
		emailError = True
	if (mainLibraryFolder==""):
		error = True
		mainLibraryFolderError = True
	if (subLibraryIDs==""):
		error = True
		subLibraryIDsError = True
	if (referenceSequences==""):
		error = True
		referenceSequencesError = True

	if not error:
		#Path to Pipeline.py
		pathToPipeline = 'secret'

		#Log file
		timeStamp = int(time.time())
		logFile =pathToPipeline+'/website/logs/logOfPipeline_'+str(timeStamp)+'.txt'

		#Pipeline command
		command = 'python3'+" "+pathToPipeline+'/pipeline.py'+' -m '+mainLibraryFolder+' -s '+subLibraryIDs+' -r '+referenceSequences+' -e '+email+' -l '+logFile
		commandString = shlex.split(command)

		#Run pipeline.py
		p = subprocess.Popen(commandString, stdout=subprocess.PIPE)

		return render_template('index.html', submissionSuccess=True, email=email, mainLibraryFolder=mainLibraryFolder, subLibraryIDs=subLibraryIDs.split(" "), referenceSequences=referenceSequences, command=" ".join(commandString))
	else:
		return render_template('index.html', submissionFailure=True, emailError=emailError, mainLibraryFolderError=mainLibraryFolderError, subLibraryIDsError=subLibraryIDsError, referenceSequencesError=referenceSequencesError)




#
#Shutdown server
#
def shutdown_server():
	func = request.environ.get('werkzeug.server.shutdown')
	if func is None:
		raise RuntimeError('Not running with the Werkzeug Server')
	func()

@app.route('/shutdown')
def shutdown():
	print("Shutting down server...")
	shutdown_server()
	return 'Server shutting down...'
#
# Pages
#
@app.route('/')
def index():
	return render_template('index.html')

@app.route('/about')
def about():
	return render_template('about.html')

if __name__ == '__main__':
	app.run(debug=True,
			port=83)
	if not app.debug:
		import logging
		from loggin.handlers import FileHandler
		file_handler = FileHandler('log.txt')
		file_handler.setLevel(logging.WARNING)
		app.logger.addHandler(file_handler)
