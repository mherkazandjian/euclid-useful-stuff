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
import migration
from migration import run_command

# some project paths
#project_path = 'EC/SGS/OU/NIR/NIRPipe/WP2500_CR'
#project_path = 'EC/SGS/SDC/I/Projects/CTEigen'
project_path = 'EC/SGS/SDC/CH/Projects/Elements'

# path to script that gets authors list. This can be obtained from
# https://bitbucket.org/atlassian/svn-migration-scripts/downloads/svn-migration-scripts.jar
path_to_svn_migration_script = '~/bin/svn-migration-scripts.jar'

# path to the svn and git repos
euclid_svn = 'http://euclid.esac.esa.int/svn'
euclid_git = 'http://euclid-git.roe.ac.uk/mkazandj'

# get the authros list
migration.get_authors(path_to_svn_migration_script, euclid_svn)

# get the project name from the project path
project_name = project_path.split('/')[-1]

# clone the repository
cmd = ('git svn clone --stdlayout --authors-file=authors.txt '
       '{euclid_svn}/{project_path} {project_name}').format(
           euclid_svn=euclid_svn,
           project_path=project_path,
           project_name=project_name)
run_command(cmd)

# add the authors file
cmd = 'git config svn.authorsfile ../authors.txt'
run_command(cmd, cwd=project_name)

# fetch the repo
cmd = 'git svn fetch'
run_command(cmd, cwd=project_name)

# clean the repo (branches/tags/remotes)
migration.create_tags(project_name)
migration.create_branches(project_name)
migration.delete_all_remotes(project_name)

# add the remotes and sync to gitlab
migration.add_remotes(project_name, euclid_git)
migration.sync(project_name)

print('done')
