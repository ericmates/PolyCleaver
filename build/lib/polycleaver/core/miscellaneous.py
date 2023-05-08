import sys

def load_bar(index, total):
    percent = 100.*(index+1)/total
    sys.stdout.write('\r')
    sys.stdout.write("Sanitizing slabs: [{:{}}] {:>3}%"
                        .format('='*int(percent/(100.0/30)),
                                30, int(percent)))