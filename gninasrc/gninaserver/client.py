#!/usr/local/bin/python

#A test script for the server, it is assumed that the receptor is pdbqt and the ligands are smina

import sys,argparse,socket,time

def recv_all(s):
    '''I can't believe that I can't find a higher level interface
    to sockets that let's me do stuff like this - read everything'''
    data = s.recv(8096)
    result = ''
    while data:
        result += data
        data = s.recv(8096)
    return result

parser = argparse.ArgumentParser(add_help=False)
parser.add_argument('-r','--receptor', help="receptor")
parser.add_argument('-l','--ligands',help="ligands")
parser.add_argument('-h','--host',help="host name")
parser.add_argument('-p','--port',help="server port")
parser.add_argument('-q','--qid',help="qid for fetch results",default=0)
parser.add_argument('-o','--out',help="out file for fetch results",default="min.sdf.gz")
args = parser.parse_args()

if args.qid > 0:
    #fetch sdf.gz for already comptued query
    s = socket.create_connection((args.host, args.port))
    s.sendall("getmols\n")
    s.sendall("%s\n" % args.qid)
    #specify filter
    s.sendall("9999999 9999999 0 0 0 0 0\n") #note that whitespace at the end is required
    
    result = recv_all(s)
    out = open(args.out,'w')
    out.write(result)
    
else:
    #open files
    rec = open(args.receptor).read()
    ligs = open(args.ligands).read()
    
    #connect to server
    s = socket.create_connection((args.host, args.port))
    
    #construct query
    s.sendall("startmin\n0\n") #cmd and oldqid
    s.sendall("receptor %d\n"%len(rec))
    s.sendall(rec)
    s.sendall("0\n") #no reorient
    s.sendall(ligs)
    
    qid = int(s.recv(1024))
    
    print qid
    #one command per a connection
    s.close()
    
    while True:
        s = socket.create_connection((args.host, args.port))
        s.sendall("getscores\n")
        s.sendall("%d\n" % qid)
        #specify filter
        s.sendall("9999999 9999999 0 0 0 0 0\n") #note that whitespace at the end is required
        
        status =  recv_all(s)
        sys.stdout.write(status)
        if  status.startswith('1'):
            break
        s.close()
        time.sleep(1)
