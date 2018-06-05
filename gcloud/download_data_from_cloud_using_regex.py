def read_file_lines(path):
    f = open(path)
    arr = f.readlines()
    for i in range(len(arr)):
        arr[i]=arr[i].rstrip()
    f.close()
    return arr

import os,re
OUT = "/Users/David/Desktop/MoTrPAC/data/pilot/"
buckets = os.system("gsutil ls > bucket_names.txt")
ignore_buckets = ["gs://motrpac-portal-projects/"]
buckets = read_file_lines('bucket_names.txt')
for b in buckets:
    if b in ignore_buckets:continue
    dirname = b.split("gs://")[1]
    os.system("mkdir "+OUT+dirname)
    os.system("gsutil cp -r " + b + " " + OUT + dirname)
