"""
<keywords>
python, migration, git, svn, main, script
</keywords>
<description>
main script that migrates a euclid repo from lousy svn to git
</description>
<seealso>
</seealso>
"""

import argparse
import re
import migration
from migration import run_command

def valid_url(url):
    """ Check if the given URL already exists
    """
    regex = re.compile(
        r'^(?:http|ftp)s?://' # http:// or https://
        r'(?:(?:[A-Z0-9](?:[A-Z0-9-]{0,61}[A-Z0-9])?\.)+(?:[A-Z]{2,6}\.?|[A-Z0-9-]{2,}\.?)|' #domain...
        r'localhost|' #localhost...
        r'\d{1,3}\.\d{1,3}\.\d{1,3}\.\d{1,3})' # ...or ip
        r'(?::\d+)?' # optional port
        r'(?:/?|[/?]\S+)$', re.IGNORECASE)

    if regex.match(url):
        return url
    else:
        raise IOError("URL '%s' is not well formatted" % url)


def migrate_project(base_svn_url, base_git_url, relative_project_url):
    """ Migrate a project from SVN to Git

    :param base_svn_url: The base url of the svn repo
    :param base_git_url: The base url of the git repo
    :param relative_project_url: The relative path of the project on the svn
     repo.
    """
    # migrate the project
    migration.migrate_project(base_svn_url=base_svn_url,
                              base_git_url=base_git_url,
                              relative_project_url=relative_project_url)


if __name__ == "__main__":
    # Create argument parser
    parser = argparse.ArgumentParser(description="Euclid migration script from SVN to Git")

    # Add arguments to parser
    parser.add_argument("base_svn_url",
                        action="store",
                        type=valid_url,
                        metavar="URL",
                        help='Base SVN URL, must be equal to "http://euclid.esac.esa.int/svn"')
    parser.add_argument("base_git_url",
                        action="store",
                        type=valid_url,
                        metavar="URL",
                        help='Base Git URL, it must begins with "http://euclid-git.roe.ac.uk"')
    parser.add_argument("relative_project_url",
                        action="store",
                        type=str,
                        metavar="PATH",
                        help='Relative project URL in the SVN repository. Eg: "EC/SGS/ST/4-2-09-TOOLS/Astromatic/Swarp"')

    # Parse command line
    args = parser.parse_args()

    # Run main function for project migration
    migrate_project(args.base_svn_url, args.base_git_url, args.relative_project_url)

