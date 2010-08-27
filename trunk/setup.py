#!/usr/bin/env python

import os

NAME = 'pylibdeconv'
AUTHOR = 'Pearu Peterson'
AUTHOR_EMAIL = 'pearu.peterson@gmail.com'
LICENSE = 'LGPL'
URL = 'http://code.google.com/p/pylibdeconv/'
DOWNLOAD_URL = 'http://code.google.com/p/pylibdeconv/downloads/list'
DESCRIPTION = 'PyLibDeconv - a Python wrapper of the Deconv library.'
LONG_DESCRIPTION = '''\
See http://code.google.com/p/pylibdeconv/ for more information.
'''
CLASSIFIERS = ''
PLATFORMS = ['Linux', 'Windows']
MAJOR               = 1
MINOR               = 0
MICRO               = 0
ISRELEASED          = not True
VERSION             = '%d.%d.%d' % (MAJOR, MINOR, MICRO)

if os.path.exists('MANIFEST'): os.remove('MANIFEST')

def write_version_py(filename='libdeconv/version.py'):
    cnt = """
# THIS FILE IS GENERATED FROM pylibdeconv/setup.py
short_version='%(version)s'
version='%(version)s'
release=%(isrelease)s

if not release:
    version += '.dev'
    import os
    svn_version_file = os.path.join(os.path.dirname(__file__),
                                   '__svn_version__.py')
    svn_entries_file = os.path.join(os.path.dirname(__file__),'.svn',
                                   'entries')
    if os.path.isfile(svn_version_file):
        import imp
        svn = imp.load_module('libdeconv.__svn_version__',
                              open(svn_version_file),
                              svn_version_file,
                              ('.py','U',1))
        version += svn.version
    elif os.path.isfile(svn_entries_file):
        import subprocess
        try:
            svn_version = subprocess.Popen(["svnversion", os.path.dirname (__file__)], stdout=subprocess.PIPE).communicate()[0]
        except:
            pass
        else:
            version += svn_version.strip()

print version
"""
    a = open(filename, 'w')
    try:
        a.write(cnt % {'version': VERSION, 'isrelease': str(ISRELEASED)})
    finally:
        a.close()


def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration(None,parent_package,top_path)
    config.add_subpackage('libdeconv')
    config.get_version('libdeconv/version.py')
    config.add_data_files(('libdeconv','LICENSE.txt'))
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup

    # Rewrite the version file everytime
    if os.path.exists('libdeconv/version.py'): os.remove('libdeconv/version.py')
    write_version_py()

    setup(
        name = NAME,
        author = AUTHOR,
        author_email = AUTHOR_EMAIL,
        license = LICENSE,
        url = URL,
        download_url = DOWNLOAD_URL,
        description = DESCRIPTION,
        long_description = LONG_DESCRIPTION,
        classifiers = filter(None, CLASSIFIERS.split('\n')),
        platforms = PLATFORMS,
        configuration=configuration)
