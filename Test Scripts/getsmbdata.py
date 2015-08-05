import urllib.request 
from smb.SMBHandler import SMBHandler

#SMB credentials
username = ""
password = "secret"

#Open URL
director = urllib.request.build_opener(SMBHandler)
fh = director.open('smb://'+username+':'+password+'@secretserver.jbei.org/miseq/secretfile')

a = fh.read()

print(a)

# Process fh like a file-like object and then close it.
fh.close()
 #No newline at end of file
