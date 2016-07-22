"""
<keywords>
python, migration, git, svn, utilities
</keywords>
<description>
module that provides implementation of function that help migrate a repo from
svn to git
</description>
<seealso>
</seealso>
"""
import os
from subprocess import Popen, PIPE
import shlex
import pdb

def run_command(cmd, *args, **kwargs):
    """wrapped around Popen

    :param str cmd: The command to be executed
    :param args: args passed to Popen
    :param kwargs: kwargs passed to Popen
    """
    print('>>> {}'.format(cmd))
    Popen(shlex.split(cmd), *args, **kwargs).wait()


def get_authors(path_to_svn_migration_script, repo_url):
    """
    dump the authors list into a file called authors.txt

    :param str path_to_svn_migration_script: path to the jar file of script
     that can be found on # https://bitbucket.org/atlassian/svn-migration-scripts/downloads/svn-migration-scripts.jar
    :param str repo_url: the url to svn repository
    """
    cmd = 'java -jar {} authors {}'.format(
        os.path.expanduser(path_to_svn_migration_script),
        repo_url)

    prcs = Popen(shlex.split(cmd), stdout=PIPE)
    stdout, stderr = prcs.communicate()

    with open('authors.txt', 'w') as fobj:
        fobj.write(stdout.decode())


def get_all_branches(project):
    """execute the command 'git branch -a' and returns the branches as a list
     of strings. The command is executed in the dir specified by the project
     parameter.

     :param str project: the name of the project
     :return: a list of strings of the branch names
    """
    cmd = 'git branch -r'
    prcs = Popen(shlex.split(cmd), cwd=project, stdout=PIPE)
    retval = []
    output, errors = prcs.communicate()

    for branch in output.encode().strip().split('\n'):
        branch_name = branch.replace('*', '').strip()
        retval.append(branch_name)
    return retval


def create_tag_from_branch(project, branch):
    """given a branch name, a tag is created and the branch is deleted. The
     command is executed in the dir specified by the project parameter.

    :param str project: the name of the project
    :param str branch: """
    run_command('git checkout {}'.format(branch), cwd=project)
    run_command('git tag {}'.format(branch.split('/')[-1]), cwd=project)


def create_branch_from_remote(project, branch):
    """given a remote branch name a local branch is created.  The command is
     executed in the dir specified by the project parameter.

    :param str project: the name of the project
    :param str branch:
    """
    run_command('git checkout {}'.format(branch), cwd=project)
    run_command('git checkout -B {}'.format(branch.split('/')[-1]), cwd=project)

def create_tags(project):
    """given a string list of branches, all branches that match pattens
     defined in this function are tagged and the branch itself is deleted
     The command is executed in the dir specified by the project parameter.
    :param str project: the name of the project
    """

    svn_tag_patterns = ['origin/tags',
                        'remotes/origin/tags']

    all_branches = get_all_branches(project)

    def is_a_tag_branch(branch):
        for pattern in svn_tag_patterns:
            if pattern in branch:
                return True
        return False
        
    tag_branches = list(filter(is_a_tag_branch, all_branches))

    for branch in tag_branches:
        create_tag_from_branch(project, branch)


def create_branches(project):
    """create branches from non tags. The command is executed in the dir
     specified by the project parameter.

    :param str project: the name of the project
    """
    svn_tag_patterns = ['origin/tags',
                        'remotes/origin/tags']

    all_branches = get_all_branches(project)

    def is_not_a_tag_branch(branch):
        for pattern in svn_tag_patterns:
            if pattern in branch:
                return False
        return True
        
    remote_branches = list(filter(is_not_a_tag_branch, all_branches))

    for remote_branch in remote_branches:
        if 'trunk' in remote_branch:
            continue
        create_branch_from_remote(project, remote_branch)


def delete_all_remotes(project):
    """deleted the remotes dir in .git/refs via "rm -fvr .git/refs/remotes/*"
    The command is executed in the dir specified by the project parameter.

    :param str project: the name of the project
    """
    run_command("rm -fvr .git/refs/remotes/origin", cwd=project)


def add_remotes(project, url_base):
    """add the remotes to the project
    git remote add origin https://<user>@bitbucket.org/<user>/<repo>.git
    The command is executed in the dir specified by the project parameter.

    :param str project: the name of the project
    :param str url_base:
    """
    run_command("git remote add gitlab {}/{}.git".format(url_base, project),
                cwd=project)

    
def sync(project):
    """
    git push -u origin --all project
    git push --tags project
    The command is executed in the dir specified by the project parameter.

    :param str project: the name of the project
    """
    run_command("git push -u gitlab --all", cwd=project)
    run_command("git push -u gitlab --tags", cwd=project)

