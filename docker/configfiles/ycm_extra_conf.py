import ycm_core
from os import getcwd
from os.path import abspath, join, isabs, normpath, exists, splitext, \
        dirname

####
# Global lists for the flags and file detection
####

##
# This is the list of default flags.
##
default_flags = [
    "-Wall",
    "-Wextra",
]

##
# C header extensions
##
c_header_extensions = [
    ".h",
]

##
# C source extensions
##
c_source_extensions = [
    ".c",
]

##
# C additional flags
##
c_additional_flags = [
    # Tell clang that this is a C file.
    "-x",
    "c",

    # Use the latest standard if possible.
    "-std=c11",
]

##
# CPP header extensions
##
cpp_header_extensions = [
    ".h",
    ".hh",
    ".H",
    ".hp",
    ".hpp",
    ".HPP",
    ".hxx",
    ".h++",
]

##
# CPP source extensions
##
cpp_source_extensions = [
    ".cp",
    ".cpp",
    ".CPP",
    ".cc",
    ".C",
    ".cxx",
    ".c++",
]

##
# CPP additional flags
##
cpp_additional_flags = [
    # Tell clang that this file is a CPP file.
    "-x",
    "c++",

    # Use the latest standard if possible.
    "-std=c++11",
]


####
# Helper functions
####

##
# Methods for file system interaction
##

def find_file_recursively(file_name, start_dir = getcwd(), stop_dir = None):
    """
    This method will walk trough the directory tree upwards
    starting at the given directory searching for a file with
    the given name.

    :param file_name: The name of the file of interest. Make sure
                      it does not contain any path information.
    :type file_name: str
    :param start_dir: The directory where the search should start.
                      If it is omitted, the cwd is used.
    :type start_dir: str
    :param stop_dir: The directory where the search should stop. If
                     this is omitted, it will stop at the root directory.
    :type stop_dir: str
    :rtype: str
    :return: The file path where the file was first found.
    """
    cur_dir = abspath(start_dir) if not isabs(start_dir) else start_dir

    while True:
        if exists(join(cur_dir, file_name)):
            # The file of interest exists in the current directory
            # so return it.
            return join(cur_dir, file_name)

        # The file was not found yet so try in the parent directory.
        parent_dir = dirname(cur_dir)

        if parent_dir == cur_dir or parent_dir == stop_dir:
            # We are either at the root directory or reached the stop
            # directory.
            return None
        else:
            cur_dir = parent_dir


def file_exists(file_name, start_dir = getcwd()):
    """
    Checks whether a file with the given file name exists in any parent
    folder of the given directory.

    :param file_name: The name of the file of interest.
    :type file_name: str
    :param start_dir: The directory where to start searching. If omitted the
                      cwd is used.
    :type start_dir: str
    :rtype: bool
    :return: True if the file was found or False if not.
    """
    return find_file_recursively(file_name, start_dir) is not None


def make_path_absolute(path, base_dir=getcwd()):
    """
    Make a given path absolute using the given base directory if it is
    not already absolute.

    :param path: The path of interest.
    :type path: str
    :param base_dir: The directory which should be used to make the
                     path absolute. If it is omitted the cwd is used.
    :type base_dir: str
    :rtype: str
    :return: The absolute path.
    """
    if isabs(path):
        return path
    else:
        return join(base_dir, path)


def script_directory():
    """
    Returns the directory where the current script is located.

    :rtype: str
    :return: The directory where the current script is located.
    """
    return dirname(__file__)


##
# Methods to check for the different source file types
##

def is_header(file_path):
    """
    Checks if the given file is a header file or not.

    :param file_path: The path to the file of interest.
    :type file_path: str
    :rtype: bool
    :return: True if the file is a header or False if not.
    """
    return is_c_header(file_path) or is_cpp_header(file_path)


def is_c_header(file_path):
    """
    Checks if the given file is a C header file or not.

    :param file_path: The path to the file of interest.
    :type file_path: str
    :rtype: bool
    :return: True if the file is a C header or False if not.
    """
    (_, extension) = splitext(file_path)

    return extension in c_header_extensions


def is_cpp_header(file_path):
    """
    Checks if the given file is a CPP header file or not.

    :param file_path: The path to the file of interest.
    :type file_path: str
    :rtype: bool
    :return: True if the file is a CPP header or False if not.
    """
    (_, extension) = splitext(file_path)

    return extension in cpp_header_extensions


def is_source(file_path):
    """
    Checks if the given file is a source file or not.

    :param file_path: The path to the file of interest.
    :type file_path: str
    :rtype: bool
    :return: True if the file is a source file or False if not.
    """
    return is_c_source(file_path) or is_cpp_source(file_path)


def is_c_source(file_path):
    """
    Checks if the given file is a C source file or not.

    :param file_path: The path to the file of interest.
    :type file_path: str
    :rtype: bool
    :return: True if the file is a C source file or False if not.
    """
    (_, extension) = splitext(file_path)

    return extension in c_source_extensions


def is_cpp_source(file_path):
    """
    Checks if the given file is a CPP source file or not.

    :param file_path: The path to the file of interest.
    :type file_path: str
    :rtype: bool
    :return: True if the file is a CPP source file or False if not.
    """
    (_, extension) = splitext(file_path)

    return extension in cpp_source_extensions


def is_c_file(file_path):
    """
    Checks if the given file is a C file or not.

    :param file_path: The path to the file of interest.
    :type file_path: str
    :rtype: bool
    :return: True if the file is a C file or False if not.
    """
    return is_c_source(file_path) or is_c_header(file_path)


def is_cpp_file(file_path):
    """
    Checks if the given file is a CPP file or not.

    :param file_path: The path to the file of interest.
    :type file_path: str
    :rtype: bool
    :return: True if the file is a CPP file or False if not.
    """
    return is_cpp_source(file_path) or is_cpp_header(file_path)


##
# Methods to manipulate the compilation flags
##

def make_absolute_flags(flags, base_dir):
    """
    Makes all paths in the given flags which are relative absolute using
    the given base directory.

    :param flags: The list of flags which should be made absolute.
    :type flags: list[str]
    :param base_dir: The directory which should be used to make the relative
                     paths in the flags absolute.
    :type base_dir: str
    :rtype: list[str]
    :return: The list of flags with just absolute file paths.
    """
    # The list of flags which require a path as next flag.
    next_is_path = [
        "-I",
        "-isystem",
        "-iquote"
    ]

    # The list of flags which require a path as argument.
    argument_is_path = [
        "--sysroot="
    ]

    updated_flags = []
    make_absolute = False

    for flag in flags:
        updated_flag = flag

        if make_absolute:
            # Assume that the flag is a path.
            updated_flag = make_path_absolute(flag, base_dir)

            make_absolute = False

        # Check for flags which expect a path as next flag.
        if flag in next_is_path:
            # The flag following this one must be a path which may needs
            # to be made absolute.
            make_absolute = True

        # Check the flags which normally expect as the next flag a path,
        # but which are written in one string.
        for f in next_is_path:
            if flag.startswith(f):
                path = flag[len(f):].lstrip()

                # Split the flag up in two separate ones. One with the actual
                # flag and one with the path.
                updated_flags.append(f)
                updated_flag = make_path_absolute(path, base_dir)

                break

        # Check for flags which expect a path as argument.
        for f in argument_is_path:
            if flag.startswith(f):
                path = flag[len(f):].lstrip()
                updated_flag = f + make_path_absolute(path, base_dir)

                break

        updated_flags.append(updated_flag)

    return updated_flags


def strip_flags(flags):
    """
    Remove leading and trailing spaces from the list of flags.

    :param flags: The list of flags which should be stripped.
    :type flags: list[str]
    :rtype: list[str]
    :return: The list of flags with leading and trailing spaces removed.
    """
    return [flag.strip() for flag in flags]


def make_final_flags(file_name, flags, base_dir = getcwd()):
    """
    Finalize the given flags for the file of interest. This step
    includes stripping the flags, making them absolute to the given
    base directory and adding the corresponding file type infos to them
    if necessary.

    :param file_name: The name of the file of interest.
    :type file_name: str
    :param flags: The flags which have been collected so far for the file.
    :type flags: list[str]
    :param base_dir: The directory which should be used to make the flags
                     absolute. If this is omitted the cwd is used.
    :type base_dir: str
    :rtype: dict[str,object]
    :return: The finalized flags for the file in the format wanted by YCM.
    """
    stripped = strip_flags(flags)
    absolute = make_absolute_flags(stripped, base_dir)

    if is_cpp_file(file_name):
        final = save_add_flags(absolute, cpp_additional_flags)
    elif is_c_file(file_name):
        final = save_add_flags(absolute, c_additional_flags)
    else:
        final = absolute

    return create_result(final)


def save_add_flags(old_flags, additional_flags):
    """
    Add additional compilation flags to an already existing list of
    compilation flags in a way that no duplication is occurring and no
    conflicting flags are added.

    As the flags which can be passed to clang are not trivial not all cases
    can be catch. However, to simplify things flags expecting an argument
    should either have the argument as next flag or separated by a space. So
    the following examples are valid:

        - "-x", "c++"
        - "-x c++"

    This is not valid and will lead to wrong behavior:

        "-xc++"

    Flags expecting an argument separated by a "=" sign should have them
    directly after the sign. So this is the correct format:

        "-std=c++11"

    :param old_flags: The list of compilation flags which should be extended.
    :type old_flags: list[str]
    :param additional_flags: The list of compilation flags which should be
                             added to the other list.
    :type additional_flags: list[str]
    :rtype: list[str]
    :return: The properly merged result list.
    """
    skip_next_af = False

    for j in range(len(additional_flags)):
        af = additional_flags[j].strip()

        argument_type = "none"
        to_add = True

        if skip_next_af:
            # The current flag is an argument for the previous flag. This
            # should have been added already.
            skip_next_af = False
            continue

        # First check if the flag has an argument as next flag.
        if len(additional_flags) > j + 1 and \
                not additional_flags[j+1].startswith("-"):
            # There is a flag after the current one which does not start with
            # "-". So assume that this is an argument for the current flag.
            argument_type = "next"
            skip_next_af = True

            af_arg = additional_flags[j+1].strip()

        # Next check if the flag has an argument separated by a " ".
        elif af.find(" ") != -1:
            # The argument for the current flag is separated by a space
            # character.
            pos = af.find(" ")

            argument_type = "next"
            af_arg = af[pos+1:].strip()
            af = af[:pos]

        # Next check if the flag has an argument separated by a "=".
        elif af.find("=") != -1:
            # The argument for the current flag is separated by a equal
            # sign.
            pos = af.find("=")

            argument_type = "same"
            af_arg = af[pos+1:].strip()
            af = af[:pos]

        # Check against all flags which are in the already contained in the
        # list.
        skip_next_of = False

        for i in range(len(old_flags)):
            of = old_flags[i].strip()

            if skip_next_of:
                # The current flag is an argument for the previous one. Skip
                # it.
                skip_next_of = False
                continue

            # If there is no argument for this flag, check for simple
            # equality.
            if argument_type == "none":
                if of == af:
                    # The flag is already in the list. So do not add it.
                    to_add = False
                    break

            # The flag is with argument. Check for these cases.

            elif argument_type == "next":
                # The argument is normally given as next argument. So three
                # cases have to be checked:
                #   1. The argument is as next flag given, too.
                #   2. The argument is separated by a space.
                #   3. The argument is directly after the flag.

                # 1
                if of == af:
                    # The flags are the same. So the arguments could be
                    # different. In any case, the additional flag should not
                    # be added.
                    to_add = False
                    skip_next_of = True
                    break

                # 2
                elif of.startswith(af) and of[len(af)] == " ":
                    # The flag is the same and the argument is separated by a
                    # space. Unimportant of the argument the additional flag
                    # should not be added.
                    to_add = False
                    break

                # 3
                elif of.startswith(af):
                    # It could be the same flag with the argument given
                    # directly after the flag or a completely different one.
                    # Anyway don't add the flag to the list.
                    to_add = False
                    break

            elif argument_type == "same":
                # The argument is normally given in the same string separated
                # by an "=" sign. So three cases have to be checked.
                #   1. The argument is in the same string.
                #   2. The argument is in the next flag but the "=" sign in
                #      the current one.
                #   3. The argument is in the next flag as well as the "="
                #      sign.

                # 1 + 2 + 3
                if of.startswith(af):
                    # 1
                    if len(of) > len(af) + 1 and of[len(af)] == "=":
                        # The argument is given directly after the flag
                        # separated by the "=" sign. So don't add the flag
                        # to the list.
                        to_add = False
                        break

                    # 2
                    elif len(of) == len(af) + 1 and of[len(af)] == "=":
                        # The argument is given in the next flag but the "="
                        # is still in the current flag. So don't add the flag
                        # to the list.
                        to_add = False
                        skip_next_of = True
                        break

                    # 3
                    elif len(of) == len(af) and len(old_flags) > i + 1 \
                            and old_flags[i+1].strip().startswith("="):
                        # The argument is given in the next flag and the "="
                        # sign is also in that flag. So don't add the flag to
                        # the list.
                        to_add = False
                        skip_next_of = True
                        break


        # Add the flags if it is not yet contained in the list.
        if to_add:
            if argument_type == "none":
                old_flags.append(af)

            elif argument_type == "next":
                old_flags.extend([af, af_arg])

            elif argument_type == "same":
                old_flags.append("{}={}".format(af, af_arg))


    return old_flags


##
# Methods to create the correct return format wanted by YCM
##

def create_result(flags, do_cache = True, **kwargs):
    """
    Create the correct return value for YCM.

    :param flags: The flags for the requested file.
    :type flags: list[str]
    :param do_cache: If the result should be cached by YCM or not. If this is
                     omitted True is used.
    :type do_cache: bool
    :param kwargs: Additional arguments.
    :type kwargs: dict[str,object]
    :rtype: dict[str,object]
    :return: A dictionary in the format wanted by YCM.
    """
    ret = {"flags": flags, "do_cache": do_cache}

    return dict(ret, **kwargs)


##
# Methods to parse the different formats supported by this script
##

def parse_compile_commands(file_name, search_base = getcwd()):
    """
    Parse the clang compile database generated by cmake. This database
    is normally saved by cmake in a file called "compile_commands.json".
    As we don't want to parse it on our own, functions provided by ycm_core
    are used. The flags corresponding to the file of interest are returned.
    If no information for this file could be found in the database, the
    default flags are used.

    :param file_name: The file for which flags should be created.
    :type file_name: str
    :param search_base: The directory at which the search for the database
                        file should start. If it is omitted the cwd is used.
    :type search_base: str
    :rtype: dict[str,object]
    :returns: The flags found in the database in the format wanted by YCM.
    """
    database_path = dirname(find_file_recursively("compile_commands.json",
            search_base))

    database = ycm_core.CompilationDatabase(database_path)

    # As headers are not in the database, we have to use the corresponding
    # source file.
    if is_header(file_name):
        (name,_) = splitext(file_name)

        # Try out all C and CPP extensions for the corresponding source file.
        for ext in (c_source_extensions + cpp_source_extensions):
            alternative_name = name + ext

            if exists(alternative_name):
                compilation_info = database.GetCompilationInfoForFile(alternative_name)

                # In the database we found flags for the alternative name
                if (compilation_info.compiler_flags_):
                    return make_final_flags(file_name, compilation_info.compiler_flags_,
                            compilation_info.compiler_working_dir_)

    elif is_source(file_name):
        compilation_info = database.GetCompilationInfoForFile(file_name)

        # We found flags for the file in the database
        if (compilation_info.compiler_flags_):
            return make_final_flags(file_name, compilation_info.compiler_flags_,
                    compilation_info.compiler_working_dir_)

    # We either don't have a proper file ending or did not find any information in the
    # database. Therefor use the default flags.
    return parse_default_flags(file_name)


def parse_clang_complete(file_name, search_base = getcwd()):
    """
    Parse the configuration file for the clang complete VIM plugin.
    Therefore it looks for a ".clang_complete" file starting at the
    given directory.

    :param file_name: The file for which flags should be created.
    :type file_name: str
    :param search_base: The directory where to start with the search for
                        the configuration file. If it is omitted the cwd is
                        used.
    :type search_base: str
    :rtype: dict[str,object]
    :returns: The flags found in the file in the format wanted by YCM.
    """
    config = find_file_recursively(".clang_complete", search_base)
    config_path = dirname(config)

    with open(config, "r") as config_file:
        flags = config_file.read().splitlines()

        return make_final_flags(file_name, flags, config_path)


def parse_default_flags(file_name):
    """
    Parse and clean the default flags to use them as result for YCM.

    :param file_name: The file for which flags should be created.
    :type file_name: str
    :rtype: dict[str,object]
    :returns: The default flags in the format wanted by YCM.
    """
    return make_final_flags(file_name, default_flags, script_directory())


####
# Entry point for the YouCompleteMe plugin
####

def FlagsForFile(file_name, **kwargs):
    """
    This method is the entry point for the YCM plugin. It is called by the
    plugin to get the all necessary compiler flags to parse a specific file
    given as argument.

    :param file_name: The path to the file for which YouCompleteMe likes to do
                      auto completion.
    :type file_name: str
    :param kwargs: Additional key word arguments.
    :type kwargs: dict[str,str]
    :rtype: dict[str,object]
    :return: The compilation flags for the file in the format wanted by YCM.
    """
    # First check for a compile_commands.json file.
    search_base = dirname(file_name)

    if file_exists("compile_commands.json", search_base):
        # There exists a compile_commands.json file. Try to use this one.
        return parse_compile_commands(file_name, search_base)
    elif file_exists(".clang_complete", search_base):
        # There exists a .clang_complete file. Try to use this one.
        return parse_clang_complete(file_name, search_base)
    else:
        # No files exists. Use the default flags.
        return parse_default_flags(file_name)

