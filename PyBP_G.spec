# -*- mode: python -*-

block_cipher = None


a = Analysis(['PyBP_G.py'],
             pathex=['/Users/Alex/code/PyBP_G'],
             binaries=[],
             datas=[('AALv.txt', '.'), ('aal_labels.txt', '.')],
             hiddenimports=['Tkinter', 'FixT', 'scipy', 'matplotlib'],
             hookspath=[],
             runtime_hooks=[],
             excludes=[],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher)
pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)
exe = EXE(pyz,
          a.scripts,
          a.binaries,
          a.zipfiles,
          a.datas,
          name='PyBP_G',
          debug=False,
          strip=False,
          upx=True,
          runtime_tmpdir=None,
          console=False , icon='pybp.icns')
app = BUNDLE(exe,
             name='PyBP_G.app',
             icon='pybp.icns',
             bundle_identifier=None)
