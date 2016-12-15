"""
<keywords>
script, python, download, simulation, ousmin, nir, ounir, sc3
</keywords>
<description>
script that downloads the simulated data. usage example:

   python download_simulations-sc3.py

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
from euclid_useful_stuff.download_ousmin_data.utils import (
    find_downloadable_links,
    find_content_directory,
    download_files_in_url)


# the url that contians data for a certain release
base_url = 'https://shared.euclid.pic.es/SC3/NIP/T2/'

# the directory where the stuff will be downloaded to
download_dir = os.path.expanduser('~/euclid/data/euclid/simulated')

# read the username and the password from a plain text file (put this file
# in a read only direcotry/file where only you and read it
# e.g the file contents should look like this (without a newline at the end)
# ----------- onlyme/ousim-dcache.txt ---------
# MY_USER_NAME MY_PASSWORD
# ---------------------------------------------
username, password = open(
    os.path.expanduser('~/onlyme/ousim-dcache.txt')).read().strip().split()

# the authentication data
auth = HTTPBasicAuth(username, password)


def download_task(base_url, task_name, download_dir, auth):
    """download the content of the task, i.e everyting recursively in:

             https://shared.euclid.pic.es/SC3/NIP/T2/34080/

    :param root_path: prefix path to the dir containing the task_name
    :param task_name: the name of the task, e.g.
    :return: list of downloaded files
    """

    task_dir = os.path.join(download_dir, task_name[1:])
    if not os.path.exists(task_dir):
        os.mkdir(task_dir)

    # generate the url of the task dir  base_url/task
    task_url = urllib2.urlparse.urljoin(base_url, task_name)
    # DEBUG
    print(task_url)

    task_subdirs = find_content_directory(task_url, auth)

    for task_subdir in task_subdirs:
        subdir_url = urllib2.urlparse.urljoin(task_url, task_subdir)
        print(subdir_url)

        subdir_abs_path = os.path.join(download_dir, task_subdir[1:])
        if not os.path.exists(subdir_abs_path):
            os.mkdir(subdir_abs_path)

        download_files_in_url(subdir_url, subdir_abs_path, auth)
        print('-'*10)


def main():
    download_files_in_url(base_url, download_dir, auth)
    find_content_directory(base_url, auth)

    for task in find_content_directory(base_url, auth):
        print(task)
        download_task(base_url, task, download_dir, auth)

if __name__ == '__main__':
    main()

print('done')
