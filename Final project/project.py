import numpy as np

def phenyl_detector(string):
    str_as_list = []
    str_as_list[:0] = string
    length = len(str_as_list)
    b = str_as_list.index('1')
    a = b-1  
    if str_as_list[a] == 'C' and str_as_list[a+1] == '1' and str_as_list[a+2] == '=' and str_as_list[a+3] == 'C' and str_as_list[a+4] == 'C' and str_as_list[a+5] == '=' and str_as_list[a+6] == 'C' and str_as_list[a+7] == 'C' and str_as_list[a+8] == '=' and str_as_list[a+9] == 'C' and str_as_list[a+10] == '1':
        str_as_list.pop(a+10)
        str_as_list.pop(a+9)
        str_as_list.pop(a+8)
        str_as_list.pop(a+7)
        str_as_list.pop(a+6)
        str_as_list.pop(a+5)
        str_as_list.pop(a+4)
        str_as_list.pop(a+3)
        str_as_list.pop(a+2)
        str_as_list.pop(a+1)
        str_as_list[a] = 'Ph'
        new_string = ''.join(str_as_list)
        return new_string
    else:
        return string

def sub_detector(string):
  i = string.index('=')
  print(i)
  length = len(string)
  print(length)
  if length >= 11:
    if '1' in string:
        string = phenyl_detector(string)
        i = string.index('=')
        length = len(string)
    else:
        string
  else:
      string
  
  if i >= length/2:
    str_as_list = []
    str_as_list[:0] = string
    str_as_list.reverse()
    for a in range(len(str_as_list)):
      if str_as_list[a] == '(':
        str_as_list[a] = ')'
      elif str_as_list[a] == ')':
        str_as_list[a] = '('
      elif str_as_list[a] == 'P':
        str_as_list[a] = 'h'
      elif str_as_list[a] == 'h':
        str_as_list[a] = 'P'
    new_string = ''.join(str_as_list)
    string = new_string

  j = string.index('=')
  if j == 1:
    sub = 'right'
    return sub

  if string[j-2] == 'C':
    if string[j+2] == 'C':
      sub = 'equal'
      return sub
    elif string[j+2] == '(':
      sub = 'right'
      return sub
    elif string[j+2] == 'P':
      sub = 'right'
      return sub
  if string[j+2] == 'C':
    if string[j-2] == ')':
      sub = 'left'
      return sub
    elif string[j-2] == 'h':
      sub = 'left'
      return sub
  if string[j-2] == 'h':
    sub = 'left'
    return sub
  if string[j+2] == 'P':
    sub = 'right'
    return sub
  if string[j-2] == ')' and string[j-2] == '(':
    if string[j-3] == 'C' and string[j+3] != 'C':
      sub = 'left'
      return sub
    elif string[j-3] != 'C' and string[j-3] == 'C':
      sub = 'right'
      return sub
    elif string[j-3] == 'C' and string[j+3] == 'C':
      if string[j-4] == 'C' and string[j+4] != 'C':
        sub = 'left'
        return sub
      elif string[j-4] != 'C' and string[j+4] == 'C':
        sub = 'right'
        return sub
      elif string[j-4] == 'C' and string[j+4] == 'C':
        sub = 'equal'
        return sub

string = 'C/C=CC(C)C'

counter = 0
for a in range(len(string)):  
  if string[a] == '/':
    counter += 1
  elif string[a] == '\\':
    counter += 1
print(counter)

for a in range(len(string)-counter):
  str_as_list = []
  str_as_list[:0] = string
  print(str_as_list)
  if str_as_list[a] == '/':
    str_as_list.pop(a)
  elif str_as_list[a] == '\\':
    str_as_list.pop(a)
  new_string = ''.join(str_as_list)
  string = new_string
  print (string)