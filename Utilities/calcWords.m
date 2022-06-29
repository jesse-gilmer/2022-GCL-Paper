function [R,L,E,W,wordMem, cellMem] = calcWords(gout,member_cut)

    
    words = [];
    for q = 1:size(gout,2)
        words(q) = bi2de(gout(:,q)'>member_cut);
        wordint(q,:) = gout(:,q)'>member_cut ;
    end
    [u_words,idu] = unique(wordint,'rows');
    
    p_words = [];
    for q = 1:size(u_words,1)
        word_in = u_words(q,:);
        p_words(q) = sum(min((wordint == word_in)'))/length(words);
    end
    nbitsm= -sum(p_words .* log(p_words));
    
    xsum = sum(u_words');
    u_words(xsum==0,:) = [];
    nwords = size(u_words,1);
    
    cMax = size(gout,2);
    max_ensemble = cMax;
    
    for q = 1:size(gout,1)
        word_membership(q) = sum(u_words(:,q));
    end
    
    
    word_membership(word_membership ==0) = [];
    mean_ensemble_membership = mean(word_membership);
    
    ensemble_strength = (size(u_words,1)/max_ensemble)/mean_ensemble_membership;
    ensemble_strength_weak = (size(u_words,1)/mean_ensemble_membership);
    
    loss_cover = 1-sum(sum(gout>member_cut)>0)/cMax;
    
    ensemble_strength_rel = (1-loss_cover) * ensemble_strength;
    R = ensemble_strength_rel;
    L = loss_cover;
    E = ensemble_strength;
    W = ensemble_strength_weak;
    wordMem = word_membership;
    cellMem = u_words;
end