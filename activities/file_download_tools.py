import requests

# Grabbed the following snippet from 
# http://stackoverflow.com/questions/16694907/how-to-download-large-file-in-python-with-requests-py
def download_file(url):
    local_filename = "../data/%s" % (url.split('/')[-1])
    # NOTE the stream=True parameter
    r = requests.get(url, stream=True)
    print "Downloading %s....." % (local_filename)
    with open(local_filename, 'wb') as f:
        for chunk in r.iter_content(chunk_size=1024): 
            if chunk: # filter out keep-alive new chunks
                f.write(chunk)
                f.flush()
    return local_filename
