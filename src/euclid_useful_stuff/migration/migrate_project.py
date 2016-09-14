"""
<keywords>
python, migration, Git, SVN, main, script
</keywords>
<description>
Main script that migrates an Euclid repo from lousy SVN to Git.

After getting the script you simply need to type:

   ~> python migrate_project.py http://euclid.esac.esa.int/svn http://euclid-git.roe.ac.uk/mkazandj EC/SGS/ST/4-2-09-TOOLS/Astromatic/Swarp --authors=authors.txt

</description>
<seealso>
</seealso>
"""

import argparse
import re
import os

import migration


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


if __name__ == "__main__":
    # Create argument parser
    parser = argparse.ArgumentParser(
        description="Euclid migration script from SVN to Git")

    # Add arguments to parser
    parser.add_argument(
        "base_svn_url",
        action="store",
        type=valid_url,
        metavar="URL",
        help='Base SVN URL, must be equal to "http://euclid.esac.esa.int/svn"')

    parser.add_argument(
        "base_git_url",
        action="store",
        type=valid_url,
        metavar="URL",
        help='Base Git URL, it must begins with'
             '"http://euclid-git.roe.ac.uk/USERNAME"')

    parser.add_argument(
        "relative_project_url",
        action="store",
        type=str,
        metavar="PATH",
        help='Relative project path in the SVN repository. Eg: '
             '"EC/SGS/ST/4-2-09-TOOLS/Astromatic/Swarp"')

    parser.add_argument(
        "--authors",
        action="store",
        type=valid_file,
        metavar="FILE",
        default="authors.txt",
        help="Authors file formatted like: login = Author Name <email>"
             "DEFAULT=authors.txt")

    parser.add_argument(
        "--project-new-name",
        action="store",
        dest='project_new_name',
        type=str,
        metavar="STRING",
        default=None,
        help="The new mae of the project dir (Git project)")

    # Parse command line
    args = parser.parse_args()

    # Run main function for project migration
    migration.migrate_project(base_svn_url=args.base_svn_url,
                              base_git_url=args.base_git_url,
                              relative_project_url=args.relative_project_url,
                              authorsfile=args.authors,
                              project_new_name=args.project_new_name)
