from __future__ import print_function
import requests
from lxml import html
import urllib2
import os

def download_file(url, local_url, auth):
    """
    taken from http://stackoverflow.com/questions/16694907/how-to-download-large-file-in-python-with-requests-py
    """
    print('downloading {} to {}'.format(url, local_url))
    request = requests.get(url, stream=True, auth=auth)
    with open(local_url, 'wb') as f:
        for chunk in request.iter_content(chunk_size=1024):
            if chunk:  # filter out keep-alive new chunks
                f.write(chunk)
    return local_url


def find_downloadable_links(url, auth):
    """given an html tree, search for all the downloadable elements and return
    them as a list"""

    # authentica and get the requested url
    request = requests.get(url, auth=auth)
    tree = html.fromstring(request.content)

    hrefs = tree.xpath('.//a[@href]')

    retval = []
    for element in hrefs:
        if 'download' in element.attrib:
            retval.append(element.attrib['href'])
    return retval


def find_content_directory(base_url, auth):
    """given an html tree, search for all the directories and return them as
     a list

     :param base_url: The url of the directory to be searched
     :param auth: The authentication object
     :return: list of directory names
    """
    # authentica and get the requested url
    request = requests.get(base_url, auth=auth)
    tree = html.fromstring(request.content)

    table = tree.xpath('.//div[2]/table')[0]

    table_rows = table.xpath('//tr')

    retval = []
    for row in table_rows:
        attribs = row.xpath('.//@href')
        if len(attribs) == 1:
            retval.append(attribs[0])
    return retval


def download_files_in_url(base_url, download_dir, auth):
    """Download files that are in the base_url in to the download_dir

    :param base_url: The base url where there might be some downloadable files
    :param download_dir: The dir where the files will be downloaded
    :param auth: The authentication to be used in the request
    """

    relative_urls = find_downloadable_links(base_url, auth)
    for relative_url in relative_urls:
        abs_url = urllib2.urlparse.urljoin(base_url, relative_url)
        local_path = os.path.join(download_dir, relative_url.split('/')[-1])
        print('downloading file:\n\t{}\nto\n\t{}'.format(
            abs_url,
            local_path
        ))
        download_file(abs_url, local_path, auth)
