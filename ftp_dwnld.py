# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 13:58:26 2015
Download and save data from ftp session with http://hermes.acri.fr/index.php?class=archive
Use for merged chlorophyll products, possibly other data
NOTE *** must change directory of subdfolder (line 17)

@author: ryan.morse
"""
import os
from ftplib import FTP
os.chdir('C://Users/ryan.morse/Desktop/Iomega Drive Backup 20171012/1 RM/3 gridded data/HERMES_monthly')
#os.chdir('C:/Users/RM/Downloads')
#os.chdir('/home/ryan/Downloads')
url='ftp.hermes.acri.fr'
ftp = FTP(url, 'ftp_hermes', 'hermes%')

#ftp://ftp.hermes.acri.fr/621314976

#ftp://ftp.hermes.acri.fr/233196241
#choose FTP subfolder
ftp.cwd('/233196241/') # Change this to folder number assigned by Hermes

listing = []
ftp.retrlines("LIST", listing.append) # find all included files

# Download to specified folder
for i in  range(0, len(listing)):
    words = listing[i].split(None, 8)
    filename = words[-1].lstrip()
    local_filename = os.path.join(r"C:\Users\\ryan.morse\Desktop\Iomega Drive Backup 20171012\1 RM\3 gridded data\HERMES_monthly", filename)
    #local_filename = os.path.join(r"C:\Users\ryan.morse\Downloads", filename)
    #local_filename = os.path.join(r"/home/ryan/Downloads", filename)
    lf = open(local_filename, "wb")
    ftp.retrbinary("RETR " + filename, lf.write, 8*1024)
    lf.close()

ftp.quit()