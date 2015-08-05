import sys, os, shutil, subprocess, time
#For getting fastq.gz and references.fasta data
from MiSeqServerData import MiSeqServerData, ICEServerData
#For command line arguments parser
from argparse import ArgumentParser
#For email function
import smtplib
from email.mime.text import MIMEText



''''''''''''''''''''''''''''''
''' SYSTEM-DEFINED VARIABLES '''
''''''''''''''''''''''''''''''

#Paths to folders and tools
pathToPipeline = "path/To/thisscript" #Path where this script, pipeline.py, is actually being run from
pathToMiSeqSequenceStorage = pathToPipeline+"/"+"SMBMiSeqData" #Where MiSeq fastq.gz sequences are stored locally
#Path to Tools
pathToTools = pathToPipeline+"/tools"
#Picard, BWA, Samtools, GATK libraries
pathToPicard = pathToTools+"/Picard"
pathToBWA = pathToTools+"/BWA"
pathToSamtools = pathToTools+"/Samtools"
pathToGATK = pathToTools+"/GATK"
# MiSeqBAMGeneration Tools, Post Processing
pathToMiSeqBAMGenerationTools = pathToTools+"/MiSeqBAMGenerationTools" #Path to ReSeq Tools
pathToPostProcessing = pathToTools+"/Postprocessing"
pathToPostProcessingScripts = pathToPostProcessing+"/scripts"
os.environ["PERL5LIB"] = pathToMiSeqBAMGenerationTools #For PERL5

#Timestamp
timeStamp = int(time.time())



''''''''''''''''''
''' FUNCTIONS  '''
''''''''''''''''''

#Formatting
def log(file, msg):
	file.write(msg+"\n")
	print(msg+"\n")
commands = [] #List of commands to set up directories and files
def doCommands(commands, writeOutput, logFile):
	startTime = time.time()
	for c in commands:
		p = subprocess.Popen(c, stdout=subprocess.PIPE)
		if (writeOutput):
			for line in p.stdout: log(logFile, line.decode())
	endTime = time.time()
	log(l, str(endTime-startTime)+" seconds")
#File commands
def mkdir(folder):
	os.makedirs(folder, exist_ok=True)
def rmdir(path):
	shutil.rmtree(path)
def cd(path):
	os.chdir(path)
#Send email
def sendEmail(email, header, content):
	msg = MIMEText(content)
	# me == the sender's email address
	# you == the recipient's email address
	msg['Subject'] = str(header)
	msg['From'] = 'MiSeq Validation Pipeline'
	msg['To'] = email
	# Send the message via our own SMTP server.
	s = smtplib.SMTP('localhost')
	s.send_message(msg)
	s.quit()


########################################################################################
#
#Parse command line arguments
#

parser = ArgumentParser()
parser.add_argument("-m", "--mainLibrary", 
					dest="mainLibrary", default="", 
                    help="Name of Folder on MiSeq Machine for the Run you want analyzed, e.g. 141212_M03257_0002_000000000-AC28N")
parser.add_argument("-s", "--subLibraries",
					nargs='+',
                    dest="subLibraries", default=[],
                    help="Sample IDs of pools you want validated, e.g. WT_1 WT_2 WT_3 32_1")
parser.add_argument("-r", "--referenceSequences", 
					nargs='+',
					dest="referenceSequences", default=[], 
                    help="Name of FASTA file where reference sequences are stored")
parser.add_argument("-e", "--email", 
					dest="email", default="", 
                    help="Email address to be notified when Pipeline is completed")
parser.add_argument("-l", "--logFile", 
					dest="logFile", default=pathToPipeline+"/website/logs/logOfPipeline_"+str(timeStamp)+".txt", 
                    help="Path to log File for pipeline's output messages")
args = parser.parse_args()
#Assign command line args to local variables
mainLibrary = args.mainLibrary
subLibraries = args.subLibraries
referenceSequences = args.referenceSequences
email = args.email
logFile = args.logFile
#Create PATH
pathToMainLibrary = pathToPipeline+"/MiSeqValidationResults/"+mainLibrary #mainLibrary, where all subLibrary folders are stored
pathToReferenceFASTA = pathToMainLibrary+"/"+"ref/references.fasta" #Reference FASTA
pathToLibrariesInfo = pathToMainLibrary+"/"+"libraries.info"

########################################################################################
#
#Set up Log File
#

#Create Log file.txt
l = open(logFile, 'w')
#Log command line arguments
log(l, "Arguments passed:")
log(l, "Main Library: "+mainLibrary)
log(l, "Sub Libraries: "+str(subLibraries))
log(l, "Reference Sequences "+str(referenceSequences))
log(l, "Email "+str(email))
log(l, "")


########################################################################################
#
#Set up folders
#

#Create Project Directory
if (os.path.isdir(pathToMainLibrary)):
 	#Folder already exists for this project, delete it and start over
 	log(l, "Folder already exists for this library. Deleting...")
 	rmdir(pathToMainLibrary)
log(l, "Folder created for "+mainLibrary)
mkdir(pathToMainLibrary)
cd(pathToMainLibrary)
#Make reference directory, ref/
mkdir("ref")


########################################################################################
#
#GET and WRITE Library metadata and FASTQ.GZ files
#

#Write files from SMB Server
log(l, "Getting sublibraries' sequence data from SMB...")
MiSeqServerThreads = []
for index, subLibraryID in enumerate(subLibraries):
	# Create new threads
	MiSeqServerThreads.append(MiSeqServerData(index, mainLibrary, subLibraryID, pathToMiSeqSequenceStorage))
	#Start thread
	MiSeqServerThreads[index].start()
# Wait for all threads to complete
for t in MiSeqServerThreads:
	t.join()
#Write metadata to libraries.info
log(l, "Writing libraries.info file...")
f = open(pathToLibrariesInfo, "w")
for t in MiSeqServerThreads:
	f.write(t.metadata)
f.close()


########################################################################################
#
#GET and WRITE Reference sequences
#

# NEEDS TO BE REWRITTEN:
# This script is given an ICE Entry's ID
# Then it queries ICE, and reads that entry's sequence
# Then it writes that sequence to a file called "references.fasta" inside the ref/ folder of the mainLibrary, as stored in the variable - "pathToReferenceFASTA"
# That's all it does, write a .fasta file to the path located in pathToReferenceFASTA

#Get files from ICE
log(l, "Getting reference sequence data from ICE...")
ICEServerThreads = []
for index, sequenceName in enumerate(referenceSequences):
	# Create new threads
	ICEServerThreads.append(ICEServerData(index, sequenceName+".fasta"))
	#Start thread
	ICEServerThreads[index].start()
# Wait for all threads to complete
for t in ICEServerThreads:
	t.join()
#Write to reference.fasta
f = open(pathToReferenceFASTA, "w")
for t in ICEServerThreads:
	f.write(t.sequence)
f.close()


''' Resulting directory structure:
	MiSeqOutput/
		seq1_r1_001.fastq.gz
		seq1_r2_001.fastq.gz
	mainLibrary/
		libraries.info
		ref/
			references.fasta (Holds FASTA file with all reference sequences)
'''


########################################################################################
#
#Generate .bam files
#

#Run prep_ref to generate .dict, .fasta.fai
log(l, "Running prep_ref...")
commands.append(["perl", pathToMiSeqBAMGenerationTools+"/prep_ref.pl", "-index", pathToReferenceFASTA, "-picard_path", pathToPicard, "-bwa_path", pathToBWA, "-samtools_path", pathToSamtools, "-bad_to_n"])
doCommands(commands, True, l)
commands = []

''' Resulting directory structure:
		mainLibrary/
			libraries.info
			ref/
				references.fasta
				references.fasta.amb
				references.fasta.ann
				references.fasta.bak
				references.fasta.bwt
				references.fasta.pac
				references.fasta.sa
				references.fasta.fai
				references.dict
'''

#Create directories for each sublibraries
log(l, "Running beta_prep_setup_dirs...")
commands.append(["perl", pathToMiSeqBAMGenerationTools+"/beta_prep_setup_dirs.pl", "-ref_fasta", pathToReferenceFASTA, "-rna", "-config", pathToLibrariesInfo])
doCommands(commands, True, l)
commands = []

''' Resulting directory structure:
		mainLibrary/
			libraries.info
			libraries.info.bak.1
			ref/
				references.fasta
				references.fasta.amb
				references.fasta.ann
				references.fasta.bak
				references.fasta.bwt
				references.fasta.pac
				references.fasta.sa
				references.fasta.fai
				references.dict
			projID_deliveryID_libName/
				bwa_dir/
					config.yml
				fastq_dir/
					symlink to seq1_r1_001.fastq.gz
					symlink to seq1_r2_001.fastq.gz
'''

#Slice sequences
log(l, "Running beta_slice_fq...")
commands.append(["perl", pathToMiSeqBAMGenerationTools+"/beta_slice_fq.pl", "-config", pathToLibrariesInfo, "-mainlibdir", pathToMainLibrary, "-reseqbindir", pathToMiSeqBAMGenerationTools])
doCommands(commands, True, l)

''' Resulting directory structure:
		mainLibrary/
			libraries.info
			libraries.info.bak.1
			ref/
				references.fasta
				references.fasta.amb
				references.fasta.ann
				references.fasta.bak
				references.fasta.bwt
				references.fasta.pac
				references.fasta.sa
				references.fasta.fai
				references.dict
			projID_deliveryID_libName/
				bwa_dir/
					config.yml
					bam_dir/
						EMPTY
					fastq_dir/
						se-32-2_S8_L001_R1_001@1.fq.gz
						se-32-2_S8_L001_R2_001@1.fq.gz
				fastq_dir/
					symlink to seq1_r1_001.fastq.gz
 					symlink to seq1_r2_001.fastq.gz
'''

#Align sliced sequences to generate .bam, .bam.bai files
log(l, "Running beta_run_alignments...")
commands.append(["perl", pathToMiSeqBAMGenerationTools+"/beta_run_alignments.pl", "-c", pathToLibrariesInfo, "-picard_path", pathToPicard, "-bwa_path", pathToBWA, "-samtools_path", pathToSamtools, "-reseqbindir", pathToMiSeqBAMGenerationTools])
doCommands(commands, True, l)
commands = []


########################################################################################
#
#Post Processing
#

#Create config.xml
log(l, "Creating config.xml for postprocessing.sh script...")
f = open("config.xml", "w")
f.write('<?xml version="1.0" encoding="UTF-8" standalone="yes"?><SBAnalysis>')
f.write('<analysis name="'+pathToMainLibrary+'" reference="'+pathToReferenceFASTA+'" location="'+pathToPipeline+'">')
for index, subLibraryID in enumerate(subLibraries):
	poolName = subLibraryID+'_libName'
	f.write('<pool name="'+poolName+'" samples="'+poolName+'"><job name="bwa_dir" protocol="bwa_dir"/></pool>')
f.write('</analysis></SBAnalysis>')
f.close()
#Postprocessing.sh
log(l, "Running postprocessing.sh...")
commands.append([pathToPostProcessing+"/postprocessing.sh", pathToPipeline, pathToMainLibrary, pathToReferenceFASTA, pathToGATK, pathToPostProcessingScripts])
doCommands(commands, True, l)
commands = []


########################################################################################
#
#Upload to ICE
#

# NOTE: The file of interest, .igv.xml, is located at:
#  pathToIGV = pathToMainLibrary+"/results/"+".igv.xml"
# NOTE: In order to get the "call" of each SubLibrary, e.g. Incomplete, Errors, Dips, Success, do:
#  calls = [] #Maps each SubLibraryID -> call name
#  for subLibraryID in enumerate(subLibraries):
#		pathToCallSummary = pathToMainLibrary+"/"+subLibraryID+"_libName"+"/bwa_dir/call_summary.txt"
#		calls[subLibraryID] = f.open(pathToCallSummary, 'r').read()
# Ask Ernst for color coding of calls.


########################################################################################
#
#Email alert
#

#Send Email notifying user that Pipeline is finished
log(l, "Email sent to "+email)
sendEmail(email, "MiSeq Validation Pipeline Finished!", "Your sequencing results for: "+mainLibrary+" and the sample IDs "+" ".join(subLibraries)+" is finished. Go to ICE to see them.")





sys.exit()


########################################################################################
################################### END ################################################
########################################################################################



''' I'm not sure if Trimgalore is necessary, or if Joel's Scripts automatically filter reads based on quality '''


''' TRIMGALORE - REMOVE ADAPTER SEQUENCES AND LOW QUALITY READS 

pathToTrimGaloreOutput = "TrimGaloreOutput" #Folder with output from Trim Galore
pathToTrimGalore = pathToTools+"/TrimGalore"
pathToCutAdapt = "/Library/Frameworks/Python.framework/Versions/2.7/bin/cutadapt" #Cutadapt command
pathToFastQC = "/usr/local/bin/fastqc"

mkdir(pathToTrimGaloreOutput)

#Construct command
commands.append([pathToTools+"/trimgalore.pl", "--quality", "20", "--illumina", "--length", "20", "--phred33", "--path_to_cutadapt", pathToCutAdapt, "--path_to_fastqc", pathToFastQC, "--fastqc", "--output_dir", pathToTrimGaloreOutput, "--paired", pathtoinputread1, pathtoinputread2])
doCommands(commands, True, l)
commands = []


Result: There is now a trimmed file in TrimGaloreOutput/ with a filtered sequence '''


