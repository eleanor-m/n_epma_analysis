''' Functions to help with the project running smoothly. '''

from pathlib import Path

def check_calczaf_folder_exists(subfolder_specified_by_user):

    # Check if folder exists, if not, ask if user to create it

    if not subfolder_specified_by_user.is_dir():
        raise FileNotFoundError(
            f'\n"{subfolder_specified_by_user}" folder does not exist. '
            'Specify the correct path to the folder you want, or create the '
            'folder, and ensure it contains a file with valence information.'
            )

def make_folder_if_it_does_not_exist(folder):
    
    if not folder.is_dir():
        create_folder = ask_user_if_folder_should_be_created(folder)
        if create_folder:
            folder.mkdir(parents=True, exist_ok=True)
            print(f'created folder {folder}')
        else:
            msg = (f'\n\nAborted - create "{folder}" or choose correct folder ' 
                    'before trying again.\n\n')
            raise FileNotFoundError(msg)   
    else:
        pass
    
    
def ask_user_if_folder_should_be_created(folder):
    check = str(input(f'\n{folder} does not exist. Do you want it to be created?'
                        ' (Y/N): ')).lower().strip()
    try:
        if check[0] == 'y':
            print('Creating folder')
            return True
        elif check[0] == 'n':
            print('Not creating folder')
            return False
        else:
            print('Invalid Input')
            return ask_user_if_folder_should_be_created(folder)
    except Exception as error:
        print("Please enter valid inputs")
        print(error)
        return ask_user_if_folder_should_be_created(folder)