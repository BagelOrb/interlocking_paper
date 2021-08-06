import os
import re

root = 'C:\\Users\\t.kuipers\\OneDrive - Ultimaker B.V\\Documents\\PhD\\interlocking_project\\paper\\paper_git'

sources = root + '\\sources'

used_files = []

tex_files = [each for each in os.listdir(root) if each.endswith('.tex')]

regex = r"\\includegraphics[^\{]*\{([^\}]*)\}"

rename = {}
for tex_file in tex_files:
    content = ""
    with open(os.path.join(root, tex_file), 'r') as file:
        content = file.read()
        r1 = re.findall(regex, content)
        used_files += r1
        for used_file in r1:
            rename[used_file] = used_file.replace('/', '-')
            content = content.replace(used_file, rename[used_file])
    with open(os.path.join(root, tex_file), 'w') as file:
        file.write(content)

used_paths = {}
for file in used_files:
    used_paths[os.path.normpath(root + "/" + file)] = file

print(used_paths)

for subdir, dirs, files in os.walk(sources):
    for file in files:
        filename = os.path.join(subdir, file)
        normpath = os.path.normpath(filename)
        if normpath in used_paths.keys():
            # print(normpath, " >>>> ", root + '\\' + rename[used_paths[normpath]])
            os.rename(normpath, root + '/' + used_paths[normpath].replace('/', '-'))