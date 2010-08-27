import os
import glob

def configuration(parent_package='',top_path=None):
    from numpy.distutils.system_info import get_info, dict_append
    from numpy.distutils.misc_util import Configuration
    package_name = 'libdeconv'
    config = Configuration(package_name,parent_package,top_path)
    deconv_src = os.path.join(config.local_path, 'src','deconv','libdeconv')
    assert os.path.isdir(deconv_src),'Directory %r not found. Fix deconv_src in %s' % (deconv_src, __file__)

    libs_info = dict (libraries = ['deconv_shared', 'gsl'])
    blas_opt = get_info('blas_opt',notfound_action=2)
    if not blas_opt:
        raise NotFoundError,'no blas resources found'
    atlas_version = ([v[3:-3] for k,v in blas_opt.get('define_macros',[]) \
                          if k=='ATLAS_INFO']+[None])[0]
    if atlas_version:
        print ('ATLAS version: %s' % atlas_version)
    dict_append(libs_info, **blas_opt)
    fftw3_info = get_info('fftw3',notfound_action=2)
    dict_append (fftw3_info, libraries = ['fftw3f'])
    if not fftw3_info:
        raise NotFoundError,'no fftw3 resources found'
    dict_append(libs_info, **fftw3_info)

    config.add_library('deconv_shared',
                       sources = glob.glob (os.path.join (deconv_src, '*.cc')),
                       depends = glob.glob (os.path.join (deconv_src, '*.h')))
    config.add_extension('_deconv',
                         sources = ['deconv.i'],
                         include_dirs = [deconv_src],
                         extra_info = libs_info
                         )

    config.make_svn_version_py()

    return config

if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(configuration=configuration)
