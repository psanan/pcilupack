#!/usr/bin/env python
import os, sys, shutil, argparse

from argparse import Namespace

try:
  WindowsError
except NameError:
  WindowsError = None

class Installer():

  def setup(self,args):
    parser = argparse.ArgumentParser(description='Install PCILUPACK')
    parser.add_argument('--src-dir', help='Location of PCILUPACK source', default=os.environ.get('PCILUPACK_DIR',os.getcwd()))
    parser.add_argument('--arch', help='Name of this configuration', default=os.environ.get('PETSC_ARCH','build'))
    parser.add_argument('--prefix', help='The installation directory', default='')
    parser.add_argument('--ranlib', help='The ranlib utility', default=os.environ.get('RANLIB',''))
    parser.add_argument('--ar-lib-suffix', help='The archive library suffix', default=os.environ.get('AR_LIB_SUFFIX',''))
    parser.parse_args(namespace=args)
    if not args.prefix:
      print '===================================================================='
      print 'No installation directory: maybe a --prefix option was not given to ./configure'
      print '===================================================================='
      sys.exit(0)
    return

  def setupDirectories(self,args):
    self.rootDir           = args.src_dir
    self.destDir           = args.prefix
    self.arch              = args.arch
    self.copies            = []
    self.ar_lib_suffix     = args.ar_lib_suffix
    self.ranlib            = args.ranlib       
    self.archLibDir        = os.path.join(self.rootDir, self.arch, 'lib')
    self.destIncludeDir    = os.path.join(self.destDir, 'include')
    self.destLibDir        = os.path.join(self.destDir, 'lib')
    self.destConfDir       = os.path.join(self.destLibDir, 'pcilupack-conf')
    return

  def checkPrefix(self):
    if not self.destDir:
      print '********************************************************************'
      print 'No installation directory: maybe a --prefix option was not given to ./configure'
      print '********************************************************************'
      sys.exit(1)
    return

  def checkDestdir(self):
    if os.path.exists(self.destDir):
      if os.path.samefile(self.destDir, self.rootDir):
        print '********************************************************************'
        print 'Incorrect prefix usage. Specified prefix same as PCILUPACK source directory'
        print '********************************************************************'
        sys.exit(1)
      if os.path.samefile(self.destDir, os.path.join(self.rootDir,self.arch)):
        print '********************************************************************'
        print 'Incorrect prefix usage. Specified destDir same as PCILUPACK build directory'
        print '********************************************************************'
        sys.exit(1)
    return

  def mkdirp(self,dst):
    if not os.path.exists(dst):
      os.makedirs(dst)
    elif not os.path.isdir(dst):
      raise shutil.Error, 'Destination is not a directory'

  def copyfile(self, src, dst, symlinks = False, copyFunc = shutil.copy2):
    dirname = os.path.dirname(dst)
    self.mkdirp(dirname)
    copyFunc(src, dst)
    return [(src,dst)]

  def copytree(self, src, dst, symlinks = False, copyFunc = shutil.copy2):
    """Recursively copy a directory tree using copyFunc, which defaults to shutil.copy2().

    The destination directory must not already exist.
    If exception(s) occur, an shutil.Error is raised with a list of reasons.

    If the optional symlinks flag is true, symbolic links in the
    source tree result in symbolic links in the destination tree; if
    it is false, the contents of the files pointed to by symbolic
    links are copied.
    """
    copies = []
    names  = os.listdir(src)
    self.mkdirp(dst)
    errors = []
    for name in names:
      srcname = os.path.join(src, name)
      dstname = os.path.join(dst, name)
      try:
        if symlinks and os.path.islink(srcname):
          linkto = os.readlink(srcname)
          os.symlink(linkto, dstname)
        elif os.path.isdir(srcname):
          copies.extend(self.copytree(srcname, dstname, symlinks))
        else:
          copyFunc(srcname, dstname)
          copies.append((srcname, dstname))
        # XXX What about devices, sockets etc.?
      except (IOError, os.error), why:
        errors.append((srcname, dstname, str(why)))
      # catch the Error from the recursive copytree so that we can
      # continue with other files
      except shutil.Error, err:
        errors.extend((srcname,dstname,str(err.args[0])))
    try:
      shutil.copystat(src, dst)
    except OSError, e:
      if WindowsError is not None and isinstance(e, WindowsError):
        # Copying file access times may fail on Windows
        pass
      else:
        errors.extend((src, dst, str(e)))
    if errors:
      raise shutil.Error, errors
    return copies

  def createUninstaller(self):
    uninstallscript = os.path.join(self.destConfDir, 'uninstall_PCILUPACK.py')
    self.mkdirp(self.destConfDir)
    f = open(uninstallscript, 'w')
    # Could use the Python AST to do this
    f.write('#!'+sys.executable+'\n')
    f.write('import os\n')

    f.write('copies = '+repr(self.copies))
    f.write('''
for src, dst in copies:
  if os.path.lexists(dst):
    os.remove(dst)
''')
    f.close()
    os.chmod(uninstallscript,0744)
    return

  def installIncludes(self):
    self.copies.extend(self.copyfile(os.path.join(self.rootDir,'PCILUPACK.h'),os.path.join(self.destIncludeDir,'PCILUPACK.h')))
    return

  def copyLib(self, src, dst):
    '''Run ranlib on the destination library if it is an archive. Also run install_name_tool on dylib on Mac'''
    import subprocess
    # Symlinks (assumed local) are recreated at dst
    if os.path.islink(src):
      linkto = os.readlink(src)
      try:
        os.remove(dst)            # In case it already exists
      except OSError:
        pass
      os.symlink(linkto, dst)
      return
    # Do not install object files
    if not os.path.splitext(src)[1] == '.o':
      shutil.copy2(src, dst)
    if os.path.splitext(dst)[1] == '.'+self.ar_lib_suffix:
      subprocess.check_call([self.ranlib,dst])
    if os.path.splitext(dst)[1] == '.dylib' and os.path.isfile('/usr/bin/install_name_tool'):
      subprocess.check_call(['/usr/bin/install_name_tool','-id',dst,dst])
    # preserve the original timestamps - so that the .a vs .so time order is preserved
    shutil.copystat(src,dst)
    return

  def installLib(self):
    self.copies.extend(self.copytree(self.archLibDir, self.destLibDir, copyFunc = self.copyLib))
    return

  def outputDestDirDone(self):
    print '''\
====================================
Installation in %s is now complete.
Now to check if the library is working do (in current directory):
make PCILUPACK_DIR=%s PETSC_ARCH="" test
====================================\
''' % (self.destDir,self.destDir)
    return

  def runsetup(self,args):
    self.setup(args)
    self.setupDirectories(args)
    self.checkPrefix()
    self.checkDestdir()
    return

  def runcopy(self):
    print '*** Installing PCILUPACK at location:',self.destDir, ' ***'
    if not os.path.exists(self.destDir):
      try:
        os.makedirs(self.destDir)
      except:
        print '********************************************************************'
        print 'Unable to create', self.destDir, 'Perhaps you need to do "sudo make install"'
        print '********************************************************************'
        sys.exit(1)
    self.installIncludes()
    self.installLib()
    return

  def rundone(self):
    self.createUninstaller()
    self.outputDestDirDone()
    return

  def run(self,args):
    self.runsetup(args)
    self.runcopy()
    self.rundone()
    return

if __name__ == '__main__':
  Installer().run(Namespace())
