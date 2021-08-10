#!/bin/bash


install_deps(){
    apt-get update
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends -o Dpkg::Options::="--force-confnew"  build-essential cmake python3-dev mono-complete golang nodejs default-jdk npm clang-tidy-9 clang-format-10 apt-transport-https ca-certificates gnupg software-properties-common wget g++-8 golang clang jsonlint jq

    pip3 install cmakelang autopep8 pylint flake8 yamllint yamlfix yamlfmt

    npm install -g npm@latest-6
    npm install --save-dev --save-exact prettier
    npm install -g fixjson

    chmod 777 /etc/vim -R

}

install_vim_82(){
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends -o Dpkg::Options::="--force-confnew"  vim vim-gtk
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
    git clone https://github.com/preservim/nerdtree.git /etc/vim/nerdtree
    git clone https://github.com/Xuyuanp/nerdtree-git-plugin.git /etc/vim/nerdtree-git-plugin
    cd /etc/vim/bundle/YouCompleteMe && git submodule update --init --recursive && python3 install.py --clang-completer --ts-completer --java-completer --cs-completer
    cd /etc/vim/bundle/vimspector && python3 install_gadget.py --enable-c --enable-cpp #--enable-python
    git clone https://github.com/lfv89/vim-interestingwords.git /etv/vim/vim-interestingwords
}

main(){
    install_deps
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends -o Dpkg::Options::="--force-confnew"  vim vim-gtk3
    install_vim_plugins
}

main
