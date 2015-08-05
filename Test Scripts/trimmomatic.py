import subprocess

pathTo_TrimmomaticJAR = "Trimmomatic-0.33/trimmomatic.jar"
pathTo_InputFASTAQ = "32-1_S7_L001_R1_001.fastq.gz"
pathTo_OutputFASTAQ = "output.fastq.gz"
pathTo_AdaptersFASTA = "Trimmomatic-0.33/adapters/TruSeq3-SE.fa"

''' Options '''
SEorPE = "SE"
phredScore = "-phred64"
IlluminaClip = "ILLUMINACLIP:"+pathTo_AdaptersFASTA
Leading = "LEADING:3"
Trailing = "TRAILING:3"
SlidingWindow = "SLIDINGWINDOW:4:15"
MinLength = "MINLEN:36"

'''  
All Options:

SLIDINGWINDOW: Perform a sliding window trimming, cutting once the average quality within the window falls below a threshold.
LEADING: Cut bases off the start of a read, if below a threshold quality
TRAILING: Cut bases off the end of a read, if below a threshold quality
CROP: Cut the read to a specified length
HEADCROP: Cut the specified number of bases from the start of the read
MINLEN: Drop the read if it is below a specified length
ILLUMINACLIP: file:<seed mismatches>:<palindrome clip threshold>:<simple clip threshold>
	seedMismatches: specifies the maximum mismatch count which will still allow a full match to be performed
	palindromeClipThreshold: specifies how accurate the match between the two 'adapter ligated' reads must be for PE palindrome read alignment.
	simpleClipThreshold: specifies how accurate the match between any adapter etc. sequence must be against a read.
'''

command = "java -jar "+pathTo_TrimmomaticJAR+" "+SEorPE+" "+phredScore+" "+pathTo_InputFASTAQ+" "+pathTo_OutputFASTAQ+" "+IlluminaClip

print(command)

p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
(output, err) = p.communicate()
