# -*- mode: python -*-

block_cipher = None


a = Analysis(['New_PyBP.py'],
             pathex=['/Users/Alex/code/PyBP'],
             binaries=[],
             datas=[('AALv.txt', '.'), ('aal_labels.txt', '.'), ('spm.surf.gii', '.')],
             hiddenimports=['traitsui.toolkit', 'traitsui', 'scipy', 'matplotlib'],
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
          name='New_PyBP',
          debug=False,
          strip=False,
          upx=True,
          runtime_tmpdir=None,
          console=False , icon='pybp.icns')
app = BUNDLE(exe,
             name='New_PyBP.app',
             icon='pybp.icns',
             bundle_identifier=None)
