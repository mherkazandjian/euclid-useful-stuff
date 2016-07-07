"""
<keywords>
script, python, download, simulation, ousmin, nir, ounir
</keywords>
<description>
script that downloads the simulated data
</description>
<seealso>
</seealso>
"""
from __future__ import print_function
import os
import requests
from requests.auth import HTTPBasicAuth
from lxml import html
import urllib2
import pdb

# the url that contians data for a certain release
base_url = 'https://shared.euclid.pic.es/SC2/NIP_R3/'

# read the username and the password from a plain text file (put this file
# in a read only direcotry/file where only you and read it
# e.g the file contents should look like this (without a newline at the end)
#----------- onlyme/ousim-dcache.txt ---------
#MY_USER_NAME MY_PASSWORD
#---------------------------------------------
username, password = open(
    os.path.expanduser('~/onlyme/ousim-dcache.txt')).read().strip().split()

# the authentication data
auth = HTTPBasicAuth(username, password)

# authentica and get the requested url
res = requests.get(base_url, auth=auth)

# get the tasks in the topdir of the relase (assuming the tasks start with
# out_EUC
tree = html.fromstring(res.content)
tasks = filter(lambda x: 'out_EUC' in x, tree.xpath('//a[@href]/text()'))

def download_file(url, local_url):
    """
    taken from http://stackoverflow.com/questions/16694907/how-to-download-large-file-in-python-with-requests-py
    """
    r = requests.get(url, stream=True, auth=auth)
    with open(local_url, 'wb') as f:
        for chunk in r.iter_content(chunk_size=1024):
            if chunk: # filter out keep-alive new chunks
                f.write(chunk)
    return local_url

def download_task(task_name):
    """

    :param task_name:
    :return:
    """

    if not os.path.exists(task_name):
        os.mkdir(task_name)

    # generate the url of the task dir  base_url/task
    url = urllib2.urlparse.urljoin(base_url, task_name)
    # DEBUG
    print(url)

    # authentica and get the requested url
    res = requests.get(url, auth=auth)

    # find the files in the url of the task
    tree = html.fromstring(res.content)
    files = filter(lambda x: 'EUC' in x, tree.xpath('//a[@title]/text()'))
    # DEBUG
    # for f in files:
    #     print(f)
    # print(len(files))

    # generate the url to the data dir: base_url/task/data/*
    data_url = urllib2.urlparse.urljoin(base_url, task_name + '/data')

    # authentica and get the requested url
    res = requests.get(data_url, auth=auth)
    tree = html.fromstring(res.content)
    data_files = filter(lambda x: 'EUC' in x, tree.xpath('//a[@title]/text()'))
    # for f in files:
    #     print(f)
    # print(len(files))


    # download the files in the task dir
    for file in files:
        link = urllib2.urlparse.urljoin(base_url, task_name + '/' + file)
        path = os.path.join(task_name, file)
        # DEBUG
        # print(link)
        # print(path)
        download_file(link, path)

    # create the data dir and download the files in the task/data dir
    data_dir_path = os.path.join(task_name, 'data')
    if not os.path.exists(data_dir_path):
        os.mkdir(data_dir_path)

    for file in data_files:
        link = urllib2.urlparse.urljoin(base_url, task_name + '/data/' + file)
        path = os.path.join(data_dir_path, file)
        # DEBUG
        # print(link)
        # print(path)
        download_file(link, path)

for task in tasks[0:1]:
    # print(task)
    download_task(task)

print('done')
