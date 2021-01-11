import os

print("Please ENTER an appropriate comment for this git commit")
comment = input()

print(comment)

cmd = "git commit -m " + '\"'  + comment  + '\"'

os.system(cmd) 
