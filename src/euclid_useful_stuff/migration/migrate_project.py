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
project_path = 'EC/SGS/ST/4-2-09-TOOLS/Astromatic/Swarp'
# project_path = 'EC/SGS/OU/NIR/NIRPipe/WP2500_CR'
# project_path = 'EC/SGS/SDC/I/Projects/CTEigen'
# project_path = 'EC/SGS/SDC/CH/Projects/Elements'

# path to script that gets authors list. This can be obtained from
path_to_svn_migration_script = '~/bin/svn-migration-scripts.jar'

# path to the svn and git repos
euclid_svn = 'http://euclid.esac.esa.int/svn'
euclid_git = 'http://euclid-git.roe.ac.uk/mkazandj'

# get the authros list
# migration.get_authors(path_to_svn_migration_script, euclid_svn)

# migrate the project
migration.migrate_project(base_svn_url=euclid_svn,
                          base_git_url=euclid_git,
                          relative_project_url=project_path)

print('done')
