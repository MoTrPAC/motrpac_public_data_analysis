import os,re
import pickle
import datetime

date_today = str(datetime.date.today())

ignored_buckets = ['gs://motrpac-portal-projects/']
file_regex = {"pilot":"pilot"}

def read_file_lines(path):
    f = open(path)
    arr = f.readlines()
    for i in range(len(arr)):arr[i]=arr[i].rstrip()
    f.close()
    return arr

def get_file_type(path):
    arr = path.split('.')
    if len(arr)<2 or path[-1]=='/' or path[-1]==':':return None
    if arr[len(arr)-1]!="gz":return arr[-1]
    return arr[-2]

# print bucket counts in a summary table
def print_bucket_type_count_table(path,bucket2count):
    o = open(path,"w")
    o.write("bucket\ttype\tcount\n")
    for b in bucket2count:
        for f in bucket2count[b]:
            o.write(b+"\t"+f+"\t"+str(bucket2count[b][f])+"\n")
    o.close()

def get_bucket_file_type_count(bucket_to_object_names,regex = None):
    bucket2count = {}
    for b in bucket_to_object_names.keys():
        b_arr = bucket_to_object_names[b]
        if len(b_arr)==0:continue
        bucket2count[b]={}
        for f in b_arr:
            if len(f)==0:continue
            if not regex is None and not re.search(regex,f,re.IGNORECASE):continue
            ftype = get_file_type(f)
            if not ftype: continue
            bucket2count[b][ftype] = bucket2count[b].get(ftype,0)+1
    return bucket2count


os.system('gsutil ls > bucket_names.txt')
buckets = read_file_lines('bucket_names.txt')
print(buckets)
bucket_to_object_names = {}
for b in buckets:
    if b in ignored_buckets:continue
    os.system('gsutil ls -r ' + b + " >tmp_objs.txt")
    bucket_to_object_names[b] = read_file_lines('tmp_objs.txt')
    print('### > ' + b)
    for a in bucket_to_object_names[b]:print(a)
#pickle.dump(bucket_to_object_names,open("bucket_to_object_names.p","wb"))

# For each bucket count the files by their type
bucket_to_object_names = pickle.load(open("bucket_to_object_names.p","rb"))
# Get all files
path = "bucket2filetypes_all_"+date_today +".txt"
bucket2count = get_bucket_file_type_count(bucket_to_object_names)
print_bucket_type_count_table(path,bucket2count)
# Get "pilot" files
path = "bucket2filetypes_pilot_"+date_today +".txt"
bucket2count = get_bucket_file_type_count(bucket_to_object_names,"pilot")
print_bucket_type_count_table(path,bucket2count)




