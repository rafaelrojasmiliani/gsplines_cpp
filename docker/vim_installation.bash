#!/bin/bash


install_deps(){
    apt-get update
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends -o Dpkg::Options::="--force-confnew"  vim vim-gtk3 vim-nox build-essential cmake python3-dev mono-complete golang nodejs default-jdk npm clang-tidy-9 clang-format-10 apt-transport-https ca-certificates gnupg software-properties-common wget g++-8 golang clang jsonlint jq
    update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-7 700 --slave /usr/bin/g++ g++ /usr/bin/g++-7
    update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-8 800 --slave /usr/bin/g++ g++ /usr/bin/g++-8
    wget -O - https://apt.kitware.com/keys/kitware-archive-latest.asc 2>/dev/null | apt-key add -
    apt-add-repository 'deb https://apt.kitware.com/ubuntu/ bionic main'
    apt-get update
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends -o Dpkg::Options::="--force-confnew"  cmake libxml2-utils clang-format clang-tidy exuberant-ctags jsonlint
    pip3 install cmakelang autopep8 pylint flake8 prospector yamllint yamlfix yamlfmt

    npm install -g npm@latest-6
    npm install --save-dev --save-exact prettier
    npm install -g fixjson

    chmod 777 /etc/vim -R

}

install_vim_82(){
    add-apt-repository ppa:jonathonf/vim -y
    apt-get update
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends -o Dpkg::Options::="--force-confnew"  vim
}

install_vim_plugins(){
    git clone https://github.com/VundleVim/Vundle.vim.git /etc/vim/bundle/Vundle.vim
    git clone https://github.com/ycm-core/YouCompleteMe.git /etc/vim/bundle/YouCompleteMe
    git clone https://github.com/vim-latex/vim-latex.git /etc/vim/bundle/vim-latex
    git clone https://github.com/preservim/tagbar.git /etc/vim/bundle/tagbar
    git clone https://github.com/jlanzarotta/bufexplorer.git /etc/vim/bundle/bufexplorer
    git clone https://github.com/dense-analysis/ale.git /etc/vim/bundle/ale
    git clone https://github.com/aklt/vim-substitute.git /etc/vim/bundle/vim-substitute
    git clone https://github.com/SirVer/ultisnips.git /etc/vim/bundle/ultisnips
    git clone https://github.com/honza/vim-snippets.git /etc/vim/bundle/vim-snippets
    git clone https://github.com/tpope/vim-fugitive.git /etc/vim/bundle/vim-fugitive
    git clone https://github.com/sukima/xmledit.git /etc/vim/bundle/xmledit
    git clone https://github.com/puremourning/vimspector.git /etc/vim/bundle/vimspector
    cd /etc/vim/bundle/YouCompleteMe && git submodule update --init --recursive && python3 install.py --clang-completer --ts-completer --java-completer --cs-completer
    cd /etc/vim/bundle/vimspector && python3 install_gadget.py --enable-c --enable-cpp #--enable-python
}

main(){
    install_deps
    install_vim_82
    install_vim_plugins
}

main
