# Routine to store methods for in-and-out processing, like reading and
# writing files

import os

import copy

########################################################################
#                              Path tools                              #
########################################################################

# Defines a function to verify if a path exists or not. If not, create
# it

def verify_path(parent_path, file_name):

    # Checks if the parent path exists

    if not os.path.exists(parent_path):

        # Initializes a list of directories to create

        directories = [""]

        # Iterates through the path

        for i in range(len(parent_path)):

            # If the character is a bar, this means the last directory
            # name has been finished

            if parent_path[-(i+1)]=="/":

                # Verifies if the last directory is not empty

                if directories[-1]!="":

                    directories.append("")

            # Otherwise, saves the characters

            else:

                directories[-1] = parent_path[-(i+1)]+directories[-1]

        # Checks if the last saved directory is empty

        if directories[-1]=="":

            directories = directories[0:-1]

        # Initializes the path

        path = ""

        # Iterates through the directories

        for i in range(len(directories)):

            # Appends this bit of directory to the path

            path += "//"+directories[-(i+1)]

            # Checks if this path exists

            if not os.path.exists(path):

                # Creates the directory

                os.mkdir(path)

    # Joins everything together

    return parent_path+"//"+file_name

########################################################################
#                              txt files                               #
########################################################################

# Defines a function to write a list into a txt file

def list_toTxt(saved_list, file_name, add_extension=True):

    # Converts the list syntax into a string using a recursion loop 
    # through the list elements

    saved_string = ""

    saved_string = recursion_listWriting(saved_list, saved_string)

    # Takes out the last comma

    saved_string = saved_string[0:-1]

    # Saves the string into a txt file

    if add_extension:

        file_name = file_name+".txt"

    txt_file = open(file_name, "w")

    txt_file.write(saved_string)

    txt_file.close()

# Defines a function to read a list from a txt file

def txt_toList(file_name):

    # Intializes the list to be read

    read_list = []
    
    # Initializes a list of indexes to inform where to append the read
    # element

    indexes_list = []

    # Initializes a counter of elements in the current sublist

    element_counterSubList = 0

    # Reads the txt file

    saved_string = ""

    try:

        with open(file_name+".txt", "r") as infile:

            saved_string = infile.read()

    except:

        raise FileNotFoundError("The file "+file_name+".txt was not fo"+
        "und while evaluating txt_toList method in file_handling_tools"+
        ".py\n")

    # Iterates through the characters of the string, but takes the first
    # and last characters out as they are the outer brackets

    # Initializes a string to store the element that is being read

    read_element = ""

    for i in range(1,len(saved_string)-1,1):

        # Retrieves the current character

        current_character = saved_string[i]

        # If the character is a comma or a closing bracket, stops the
        # reading of the element

        if current_character=="," or current_character=="]":

            if len(read_element)>0:

                # Converts the element

                try:

                    # Tries to convert it to an integer

                    read_element = int(read_element)

                except:

                    # Tries to convert it to a float

                    try:

                        read_element = float(read_element)

                    except:

                        print("Could not convert", read_element, "to a"+
                        " number\n")

                # Appends the element to the list using the list of in-
                # dexes

                #print("\nAdds element to list")

                read_list = recursion_listAppending(read_element, 
                read_list, indexes_list, -len(indexes_list))

                #print("\nUpdated list:", read_list, "\n")

                # Clears the read element

                read_element = ""

            # If the character is a comma, updates the number o elements
            # in the current sublist

            if current_character==",":

                element_counterSubList += 1

            # If the character is a closing bracket, takes the last in-
            # dex out of the list of indexes and updates the counter of
            # elements in the sublist to the remaining last index

            if current_character=="]" and len(indexes_list)>1:

                indexes_list = indexes_list[0:-1]

                element_counterSubList = indexes_list[-1]

        # If the current character is an opening bracket, it sinalizes a
        # new sublist, hence, updates list of indexes and the counter of
        # elements in the sublist

        elif current_character=="[":

            #print("\nAdds new sublist")

            # Adds a new empty list to the read list
            
            read_list = recursion_listAppending([], read_list, 
            indexes_list, -len(indexes_list))

            #print("\nUpdated list:", read_list, "\n")

            # Updates the list of indexes and the counter of elements in
            # the sublist (makes it 0, as a new empty list is added)

            indexes_list.append(element_counterSubList)

            element_counterSubList = 0

        # Otherwise, it must be a valid element to be read

        else:

            read_element += current_character

    # Returns the read list

    return read_list

########################################################################
#                           Recursion tools                            #
########################################################################

# Defines a function to recursively access the elements of a list and
# write the list syntax into a string

def recursion_listWriting(accessed_list, saved_string):

    # Adds the initial bracket as THIS list has just begun

    saved_string += "["

    # Iterates through the element of the list

    for i in range(len(accessed_list)):

        # Accesses the i-th element

        accessed_element = accessed_list[i]

        # Tests if the element is another list. If it is, recalls the 
        # recursion function to start the process all over again

        if isinstance(accessed_element, list):

            saved_string = recursion_listWriting(accessed_element, 
            saved_string)

        # Otherwise, returns the value, as the stopping criterion of the
        # recursion is precisely when the element is no longer another 
        # list

        else:

            # Updates the string

            saved_string += str(accessed_element)+","

    # Takes out the last character, for it is the last comma and; then, 
    # adds the bracket as THIS list has ended

    if saved_string[-1]==",":

        saved_string = saved_string[0:-1]

    saved_string += "],"

    # Returns the string

    return saved_string

# Defines a function to recursively access a list using a list of inde-
# xes to append an element to the list at the last index

def recursion_listAppending(added_element, accessed_list, indexes_list, 
indexes_counter=None):
    
    #print("Added element:", added_element)

    #print("Accessed list:", accessed_list)

    #print("Indexes list:", indexes_list)

    #print("Indexes counter:", indexes_counter)

    # Verifies if it was not given

    if indexes_counter==None:

        indexes_counter = -len(indexes_list)

    # Verifies whether the index counter is negative

    elif indexes_counter>0:

        indexes_counter = -1*indexes_counter

    # The indexes counter is a counter for the list of indexes, it is 
    # negative. This works well because at each index of the accessed 
    # list, the counter is added 1, thus, when it is zero, it sinalizes
    # that the position has been found, where the new element must be 
    # added

    if indexes_counter==0:

        accessed_list.append(added_element)
    
    # Otherwise, gets the list at the desired index, sends it further to
    # the recursion and brings it back rightfully appended

    else:

        sub_list = accessed_list[indexes_list[indexes_counter]]

        # Updates the index counter

        indexes_counter += 1

        # Recovers the appended sub_list and allocates it into the ori-
        # ginal list

        accessed_list[indexes_list[indexes_counter-1]] = recursion_listAppending(
        added_element, sub_list, indexes_list, indexes_counter=
        indexes_counter)

    # Returns the corrected list

    return accessed_list

# Defines a function to get the combinations of possible indexes given
# a list of the number of possible indexes per index

def get_indexesCombinations(n_indexesList):

    return get_indexesCombinationationsRecursion(n_indexesList, 
    index_combinations=[])

# Defines a function to do the recursion activity for the construction 
# of the list of possible index combinations

def get_indexesCombinationationsRecursion(n_indexesList, 
index_combination=[], index_combinations=[], component=0):

    if len(index_combination)==0:

        index_combination = [0 for i in range(len(n_indexesList))]

    # Verifies if the component if the component is larger than the 
    # length of the list

    if component>=len(n_indexesList):

        # Adds another combination

        index_combinations.append(copy.deepcopy(index_combination))

        return index_combinations
    
    # Otherwise, iterates through the indexes

    else:

        for i in range(n_indexesList[component]):

            # Adds this index

            index_combination[component] = i+0

            # Sends the index combination down the line

            index_combinations = get_indexesCombinationationsRecursion(
            n_indexesList, index_combination=index_combination, 
            index_combinations=index_combinations, component=(component+
            1))

    return index_combinations

########################################################################
#                               Testing                                #
########################################################################

def general_test():

    t=[1,[2,[3,4], 5, [0,10]]]

    print("Original t:", t)

    indexes_counter = [1, 1]

    t = recursion_listAppending(6, t, indexes_counter)

    list_toTxt(t, "test")

    print("Appended t:", t)

    t = txt_toList("test")

    print("Read t:    ", t)

def test_indexBuilder():

    n_indexesList = [3]

    indexes = get_indexesCombinations(n_indexesList)

    print(len(indexes), "combinations:", indexes, "\n\n")

    n_indexesList = [3,3]

    indexes = get_indexesCombinations(n_indexesList)

    print(len(indexes), "combinations:", indexes, "\n\n")

    n_indexesList = [3,3,3]

    indexes = get_indexesCombinations(n_indexesList)

    print(len(indexes), "combinations:", indexes, "\n\n")

#general_test()

#test_indexBuilder()