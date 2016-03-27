#!/usr/bin/env python2
#-*- coding: utf-8 -*-

##This software (including its Debian packaging) is available to you under the terms of the GPL-3, see "/usr/share/common-licenses/GPL-3".
##Software is created and maintained by Laboratory of Biomolecular Systems Simulation at University of Gdansk.
##Contributors:
##- Tomasz Makarewicz (btchtm@ug.edu.pl)

import os, shutil, tarfile
import socket

os.chdir(os.getenv("HOME"))
project_name = 'project1'
dynamics_dir = os.getenv("HOME")+'/.dynamics/'
project_dir = dynamics_dir+project_name + '/'
##Clean "project1" temporary directory if present.
try:
	shutil.rmtree(project_dir)
except:
	pass
##Creating "project1" temporary directories
if os.path.isdir(project_dir) == False:
	os.makedirs(project_dir)

##Receive data
# Create a socket object
s = socket.socket()
# Get local machine name
host = socket.gethostname()
# Reserve a port for your service.
port = 60000

s.connect((host, port))
s.send("Hello server!")

with open(dynamics_dir + 'project1.tar.bz2', 'wb') as f:
	print 'file opened'
	while True:
		print('receiving data...')
		data = s.recv(1024)
		print('data=%s', (data))
		if not data:
			break
		# write data to a file
		f.write(data)
		f.close()

print('Successfully get the file')
s.close()
print('connection closed')

##Perform calculations
tar = tarfile.open(dynamics_dir + 'project1.tar.bz2')
tar.extractall(project_dir)
tar.close()
os.remove(dynamics_dir + 'project1.tar.bz2')
#Dynamics
tar = tarfile.open(dynamics_dir + 'project2.tar.bz2', "w:bz2")
tar.add(project_dir, recursive=True, arcname="project2")
tar.close()

##Sends results back
# Reserve a port for your service.
port = 60000
# Create a socket object
s = socket.socket()
# Get local machine name
host = socket.gethostname()
# Bind to the port
s.bind((host, port))
# Now wait for client connection.
s.listen(5)

print 'Server listening....'

while True:
	# Establish connection with client.
	conn, addr = s.accept()     
	print 'Got connection from', addr
	data = conn.recv(1024)
	print('Server received', repr(data))

	filename='mytext.txt'
	f = open(filename,'rb')
	l = f.read(1024)
	while (l):
		conn.send(l)
		print('Sent ',repr(l))
		l = f.read(1024)
		f.close()

	print('Done sending')
	conn.send('Thank you for connecting')
	conn.close()
