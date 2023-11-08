

#### TKSM AUTOCOMPLETE ####

_tksm_auto_comp(){
    cur="${COMP_WORDS[COMP_CWORD]}"
    submodule="${COMP_WORDS[1]}"
    if [[ $COMP_CWORD == 1 ]]; then
        COMPREPLY=($(compgen -W "$(tksm list)" -- "${submodule}"));
    elif [[ "$cur" == --* ]]; then
        COMPREPLY=($(compgen -W "$(tksm ${submodule} --list | awk '{print "--"$0;}')" -- "${cur}"))
    else
        COMPREPLY=()
    fi
}
complete -o default  -F _tksm_auto_comp tksm
###\ TKSM AUTOCOMPLETE ###
