import csv, requests, cStringIO, re, os

USER = 'michalsa@tx.technion.ac.il'
PASSWORD = 'wallpaper947jgi'

def download(project):
    
    # Fill in your details here to be posted to the login form.
    payload = {
        'login': USER,
        'password': PASSWORD
    }
    
    chunk_size = 1024
    
    # Use 'with' to ensure the session context is closed after use.
    with requests.Session() as s:
        s.post('https://signon.jgi.doe.gov/signon/create', data=payload)
        # print the html returned or something more intelligent to see if it's a successful login page.
        #print s
    
        # An authorized request.
        r = s.get('http://genome.jgi.doe.gov/ext-api/downloads/get-directory?organism=%s' % project)
        for group in re.search('<folder name="(\w+)">', r.text).groups():
            print group
            files = re.search('<folder name="%s">\n(<file.*>\n)+</folder>' % group, r.text).group(0)
            fileTags = re.findall('(<file.*>)\n', files)
            for fileTag in fileTags: #print fileTag
                print fileTag
                if re.search(' filename="(\w+\.tar\.gz)"', fileTag):
                    fileName = re.search(' filename="(\w+\.tar\.gz)"', fileTag).group(1)
                    fileData = [] 
                    [fileData.append(el) for el in re.search('url="(.*)[&amp;(url=.*)]?\.tar\.gz"', fileTag).groups()]
                    fileUrl = fileData[0] 
                    if len(fileData) > 1:
                        fileParam = fileData[1]+'.tar.gz'
                    else:
                        fileUrl+='.tar.gz'
                        fileParam = None
                    #print fileUrl, fileParam
                    
                    #print "".join(['http://genome.jgi.doe.gov', fileUrl])
                    
                    fileReq = s.get("".join(['http://genome.jgi.doe.gov', fileUrl]), params=fileParam, stream=True)
                    print fileReq.url
                    #print fileReq.headers['content-type']
                                
                    if project+'.tar.gz' in os.listdir('/home/michalsa/workspace/python'):
                        print "file already exists in save directory"
                    
                    else:
                        if fileReq.headers['content-type'] in ['application/x-gzip', 'application/octet-stream']:
                            with open(project+'.tar.gz', 'wb') as fd:
                                #fd.write(fileReq.content)
                                
                                for chunk in fileReq.iter_content(chunk_size):
                                    #print chunk
                                    if chunk:
                                        #print chunk
                                        fd.write(chunk)
                                        fd.flush()
                                        #os.fsync(fd.fileno())            
                                
                        else:
                            print fileReq.content
                            #return False
                            #return "No IMG Data file for %s" % project[-1]
                    
                    return fileName
                
    print "no tarfile"
    return False
        
if __name__ == "__main__":
    projectsFile = '/home/michalsa/Documents/research/humanGut_genome-projects.csv'
    projectsHandle = open(projectsFile, 'r')
    projectsData = list(csv.reader(projectsHandle))
    for project in projectsData[1:]: 
        if project[0]:
            print project[0]
            status = download('IMG_%s' % project[0])
            print status
            if status == False:
                break

