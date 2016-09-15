#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Script to migrate a Euclid SVN repo to the new Git repo.

Usage:
    ./migrate_project.py --relative-project-url <path/in/svn/repo>
                         --base-git-repo <username>
                         --project-new-name <name>
                         --authors <author_file.txt>

Example:
    ./migrate_project.py --relative-project-url EC/SGS/ST/4-2-09-TOOLS/Astromatic/Swarp
                         --base-git-repo ST-TOOLS
                         --project-new-name Swarp
                         --authors author.txt

Keywords:
    python, migration, Git, SVN, script

"""

import argparse
import re
import os

import migration
from os.path import join

BASE_SVN_URL = "http://euclid.esac.esa.int/svn"
BASE_GIT_URL = "http://euclid-git.roe.ac.uk"


def valid_url(url):
    """
    Checks if the given URL is a valid URL HTTP/FTP[s]://somethig

    :param url: An URL string to be tested
    """
    # Django based URL regex
    regex = re.compile(
        r'^(?:http|ftp)s?://'  # http:// or https://
        r'(?:(?:[A-Z0-9](?:[A-Z0-9-]{0,61}[A-Z0-9])?\.)+(?:[A-Z]{2,6}\.?|'  # -
        r'[A-Z0-9-]{2,}\.?)|'  # ..>domain...
        r'localhost|'  # localhost...
        r'\d{1,3}\.\d{1,3}\.\d{1,3}\.\d{1,3})'  # ...or ip
        r'(?::\d+)?'  # optional port
        r'(?:/?|[/?]\S+)$', re.IGNORECASE)

    # Test the URL and raise an error if doesn't match the regular expression
    if regex.match(url):
        return url
    else:
        raise IOError("URL '{}' is not well formatted".format(url))


def valid_file(fname):
    """
    Checks if the given file name is valid and return the absolute path to the
    file.

    :param fname: File path to be checked
    """
    if os.path.isfile(fname):
        return os.path.abspath(fname)
    else:
        raise IOError("No such file '{}'".format(fname))


def parse_args():
    """Create argument parser"""
    parser = argparse.ArgumentParser(
        description="Euclid migration script from SVN to Git")

    parser.add_argument(
        "--base-git-repo",
        dest='base_git_repo',
        type=str,
        required=True,
        help='Base Git repository: username OR group name')

    parser.add_argument(
        "--relative-project-url",
        dest="relative_project_url",
        type=str,
        required=True,
        help='Relative project path in the SVN repository. Eg: '
             '"EC/SGS/ST/4-2-09-TOOLS/Astromatic/Swarp"')

    parser.add_argument(
        "--authors",
        type=valid_file,
        default="authors.txt",
        help="Authors file formatted like: login = Author Name <email>"
             "DEFAULT=authors.txt")

    parser.add_argument(
        "--project-new-name",
        dest='project_new_name',
        type=str,
        default=None,
        help="The new mae of the project dir (Git project)")

    return parser.parse_args()


def main():
    """Main program"""
    # Parse args
    args = parse_args()
    # Run main function for project migration
    migration.migrate_project(base_svn_url=BASE_SVN_URL,
                              base_git_url=join(BASE_GIT_URL,
                                                args.base_git_repo),
                              relative_project_url=args.relative_project_url,
                              authorsfile=args.authors,
                              project_new_name=args.project_new_name)


if __name__ == "__main__":
    main()
