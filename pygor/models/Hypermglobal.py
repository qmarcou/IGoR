class Hypermutationglobal

def read_hyper_model(filename):
    with open(filename,'r') as file:
        line = file.readline()
        while(line.find('#Hypermutationglobalerrorrate')==-1):
            line=file.readline()
            
        semi_colon_index = line.find(';')
        next_semicolon_index = line.find(';',semi_colon_index+1)
        nmer_size = int(line[semi_colon_index+1:next_semicolon_index])
        ei_contributions = zeros((nmer_size,4))
        
        line = file.readline()
        strip_line = line.rstrip('\n') #Remove end of line character
        R_val = float(strip_line)
        
        line = file.readline()
        strip_line = line.rstrip('\n') #Remove end of line character
        semi_colon_index = -1
        next_semicolon_index = strip_line.find(';')
        
        current_nucl = 0
        current_pos = 0
        
        ended = False
        
        while (next_semicolon_index!=-1):
            
            ei_contributions[current_pos,current_nucl] = float(strip_line[semi_colon_index+1:next_semicolon_index])
            
            line = file.readline()
            semi_colon_index = next_semicolon_index
            next_semicolon_index = strip_line.find(';',semi_colon_index+1)
            
            if(current_nucl<3):  
                current_nucl+=1
            else:
                current_nucl=0
                current_pos+=1
        ei_contributions[nmer_size-1,3] = float(strip_line[semi_colon_index+1:])
                
    return [nmer_size,R_val,ei_contributions]

def compute_expected_gene_mut_freq(gene_index,event_nickname,genmodel,hyperm_Nmer_len,hyperm_Pmut):
    for real in genmodel.get_event(event_nickname,by_nickname=True).realizations:
        if real.index == gene_index:
            gene_seq = real.value
            
    expected_coverage_list = []
    
    nt_seq = []
    for i in range(0,len(gene_seq)):
        
        if i>=hyperm_Nmer_len:
            #Remove first nucleotide/coverage value
            nt_seq.pop(0)
            
        nt_seq.append(nt2int(gene_seq[i]))
        
        if len(nt_seq)==hyperm_Nmer_len :
            expected_coverage_list.append(hyperm_Pmut[tuple(nt_seq)][0])
            
        else:
            if i>=(hyperm_Nmer_len-1)/2:
                expected_coverage_list.append(np.nan)
                
    expected_coverage_list += [np.nan]*((hyperm_Nmer_len-1)/2)
    
    return expected_coverage_list

#General recursive method to compute/do anything on every possible Nmers
def recurs_Nmer_apply_inner(current_pos,current_Nmer,Nmer_length,output_array,func_to_apply,*args):
    for i in range(0,4):
        new_Nmer = current_Nmer[:]
        new_Nmer.append(i)
        if(current_pos<(Nmer_length-1)):
            recurs_Nmer_apply_inner(current_pos+1,new_Nmer,Nmer_length,output_array,func_to_apply,*args)
        else:
            output_array[tuple(new_Nmer)] = func_to_apply(new_Nmer,*args)

def recurs_Nmer_apply(Nmer_length,func_to_apply,*args):
    output_array = empty(4*ones(Nmer_length))
    recurs_Nmer_apply_inner(0,[],Nmer_length,output_array,func_to_apply,*args)
    return output_array
    

def compute_P_mut_Nmers(Nmer,hyperm_model):
    if len(Nmer) == hyperm_model[0]:
        score=hyperm_model[1]
        for i,int_nt in zip(range(0,hyperm_model[0]),Nmer):
            score*=exp(hyperm_model[2][i,int_nt])
        proba = score/(1+score)
        #print(proba)
        return proba

def plot_hyperm_contributions(hyperm_model,ax,display_legend = True , display_left_label=True , display_bottom_label_ticks = True,color_set=['#79E041','#0040FF','#FBD632','#FF0000']):
    max_min_array = zeros((2,hyperm_model[0]))
    hyperm_contributions = asarray(hyperm_model[2])
    hyperm_half_len = (hyperm_model[0]+1)/2
    legend_list = []
    width = .75
    for i in range(0,4):
        col = color_set[i]
        bool_mask = hyperm_contributions[:,i]>=0
        rev_bool_mask = hyperm_contributions[:,i]<0
        rect_neg = ax.bar(find(rev_bool_mask)+1-width/2 - hyperm_half_len,hyperm_contributions[:,i][rev_bool_mask],bottom =max_min_array[1][rev_bool_mask],color=col,linewidth=0.1) #Plot the negative oens
        rect_pos = ax.bar(find(bool_mask)+1-width/2 - hyperm_half_len,hyperm_contributions[:,i][bool_mask],bottom =max_min_array[0][bool_mask],label=i,color=col,linewidth=0.1) #Plot the positive ones
        max_min_array[0][bool_mask] += hyperm_contributions[:,i][bool_mask]
        max_min_array[1][rev_bool_mask] += hyperm_contributions[:,i][rev_bool_mask]
        if any(bool_mask):
            legend_list.append(rect_pos[0])
        else:
            legend_list.append(rect_neg[0])
    if display_legend:
        ax.legend(tuple(legend_list) , ('A','C','G','T'),loc="upper left",frameon=legend_frameon,
                  fontsize = legend_fontsize , markerscale = legend_markerscale , borderaxespad = legend_markerscale,ncol=4,handlelength = .5 ,
                  handletextpad = .25 , columnspacing = .5)
    ax.hlines(0,-(hyperm_half_len),(hyperm_half_len))
    ax.set_yticks(linspace(-1.0,1.0,21),minor=True)
    ax.set_xticks(range(-hyperm_half_len+1,hyperm_half_len))
    ax.set_ylim(-.99,.99)
    ax.set_xlim((-(hyperm_half_len-.1),(hyperm_half_len-.1)))
    
    if display_left_label:
        ax.set_ylabel(r"$e_i(\sigma)$",fontdict = label_title_font,labelpad=axis_title_pad-5)
    
    if display_bottom_label_ticks:
        ax.set_xlabel("Position",fontdict = label_title_font,labelpad=axis_title_pad)
    else:
        ax.tick_params(labelbottom=False)
        
        
    ax.tick_params(axis='x',which='both',labelsize=tick_label_fontsize,top='off')
    ax.tick_params(axis='y',which='major',length = tick_len*(2.0/3.0) , width = tick_width/2.0 , labelsize = tick_label_fontsize ,right='off')
    ax.tick_params(axis='y',which='minor',length = tick_len/2.0 , width = tick_width/2.0,right='off')
        
    return ax
