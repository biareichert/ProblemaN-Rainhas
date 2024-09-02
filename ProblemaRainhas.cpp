#include <bits/stdc++.h>
#include <vector>
/* Anotações

    Autor: Geremias Corrêa, Beatriz Reichert
    Discplina: OCEV
    Semestre: 2021/1

    Exemplo de execução: g++ -Wall -O3 -o exec nomeCodigo.cpp
        ./exec <COD> <RUN> <GEN> <POP> <DIM>
*/

using namespace std;

// Realiza uma distribuição uniforme de valores entre 0 a 1 (0.00 a 1.00).
random_device randomm;
mt19937 engine{randomm()};
uniform_real_distribution<double> distribuicaoBool{0.0, 1.0};
uniform_real_distribution<double> distribuicaoInt{-5.0, 10.0};
uniform_real_distribution<double> distribuicaoReal{-10.0, 10.0};

void printPopulation(vector<vector<double>> pop, int popSize, int chromosomeSize, string message){
    cout << endl << message << endl;
    for (int i = 0; i < popSize; i++){
        cout << "Cr.["<< i << "]: [";
        for (int j = 0; j < chromosomeSize; j++){
            if (j == chromosomeSize - 1) cout << pop[i][j];
            else  cout << pop[i][j] << ", ";
        }
        cout << "]" << endl;
    }

    cout << endl;
}

void printFellow(vector<double> fellow, int popSize, int chromosomeSize, string message){
    cout << endl << message << endl;
    cout << "Cr.: [";
    for (int j = 0; j < chromosomeSize; j++){
        if (j == chromosomeSize - 1) cout << fellow[j];
        else  cout << fellow[j] << ", ";
    }
    cout << "]" << endl;
}

vector<vector<double>> gerarConfiguracaoAleatoria(int codificationType, int executionsNumber, int numberGenerations, int populationSize, int chromosomeSize){
    vector<double> vetor_auxiliar;
    vector<vector<double>> vetor_final;
    uniform_real_distribution<double> distribuicaoPerm{0.0, (double) chromosomeSize};
    int cont;

    switch (codificationType){
    case 0: //binária
        for(int i = 0; i < populationSize; i++){
            for(int j = 0; j < chromosomeSize; j++){
                double aleatorio = distribuicaoBool(engine);
                if(aleatorio > 0.5)
                    vetor_auxiliar.push_back(1);
                else
                    vetor_auxiliar.push_back(0);
            }
            vetor_final.push_back(vetor_auxiliar);
            vetor_auxiliar.clear();
        }
        break;
    case 1: //int 
        for(int i = 0; i < populationSize; i++){
            for(int j = 0; j < chromosomeSize; j++){
                int aleatorio = (int) distribuicaoInt(engine);
                vetor_auxiliar.push_back(aleatorio);
            }
            vetor_final.push_back(vetor_auxiliar);
            vetor_auxiliar.clear();
        }
        break;
    case 2: //int-perm
        for(int i = 0; i < populationSize; i++){
            for(int j = 0; j < chromosomeSize; j++){
                cont = 0;
                while(1) {
                    int aleatorio = (int) distribuicaoPerm(engine);
                    for (int k : vetor_auxiliar){
                        if (aleatorio == k){
                            cont++;
                            break;
                        }
                    }
                    if (cont == 0){
                         vetor_auxiliar.push_back(aleatorio);
                         break;
                         
                    }
                    cont = 0;
                }
            }
            vetor_final.push_back(vetor_auxiliar);
            vetor_auxiliar.clear();
        }
        break;
    case 3: //real
        for(int i = 0; i < populationSize; i++){
            for(int j = 0; j < chromosomeSize; j++){
                double aleatorio = distribuicaoReal(engine);
                vetor_auxiliar.push_back(aleatorio);
            }
            vetor_final.push_back(vetor_auxiliar);
            vetor_auxiliar.clear();
        }
        break;

    default:
        return vetor_final;
    }

    return vetor_final;
}

//link function: http://www.vision.ime.usp.br/~pmiranda/mac122_2s14/aulas/aula20/aula20.html
int howManyQueensCantBeAttacked(vector<double> linhas, int tam){
    int i;
    int x,y;
    int xx,yy;
    int cont = tam - 1, cont_final = 0;
    char vet_attack[tam];

    for (i = 0; i < tam; i++){
        vet_attack[i] = '0'; //iniciando o vetor com nenhuma rainha sendo atacada
    }

    for(i = 0; i < tam; i++){
        x = i;
        y = linhas[i];
    
        xx = x;
        yy = y;
        while(1){
            xx += 1;
            yy -= 1;
            if(xx > tam - 1 || yy < 0) break;
            
            if(yy == linhas[xx]) {
                cont -= 1;
                vet_attack[x] = 'A'; //rainha sendo atacada
                vet_attack[xx] = 'A'; //rainha sendo atacada
                break;
            }
        }
    
        xx = x; //tirar
        yy = y;
        while(1){
            xx -= 1; // +1
            yy -= 1; // + 1
            if(xx < 0 || yy < 0) break;
            
            if(yy == linhas[xx]) {
                cont -= 1;
                vet_attack[x] = 'A'; //rainha sendo atacada
                vet_attack[xx] = 'A'; //rainha sendo atacada
                break;
            }
        }
    }

    for(int i=0; i < tam; i++){
        if(vet_attack[i] == '0'){
            cont_final += 1; //contando quantas rainhas não foram atacadas
        }
    }
    return cont_final;
}

vector<vector<int>> numerar_tabuleiro(int chromosomeSize){
    vector<vector<int>> tabuleiro;
    vector<int> auxiliar;

    int cont = 1;

    for(int i = 0; i < chromosomeSize; i++){
        for(int j = 0; j < chromosomeSize; j++){
            auxiliar.push_back(cont);
            cont++;
        }
        tabuleiro.push_back(auxiliar);
        auxiliar.clear();
    }

    return tabuleiro;
}

float funcObjetivo(vector<double> individuo, int chromosomeSize){
    vector<vector<int>> tabuleiro_enumerado;

    tabuleiro_enumerado = numerar_tabuleiro(chromosomeSize);
    
    float somaValores = 0;

    for(int i = 0; i < chromosomeSize; i++){
        if(i % 2 == 0){
            somaValores = somaValores + sqrt(tabuleiro_enumerado[i][int(individuo[i])]);
        }else{
            somaValores = somaValores + log10(tabuleiro_enumerado[i][int(individuo[i])]);
        }
    }

    return somaValores;
}

float maxZToOurContext(int chromosomeSize){
    int greaterElementFromColumn = chromosomeSize;
    float maxZ = 0;

    for (int i = 0; i < chromosomeSize; i++){
        if (i % 2 == 0){
            maxZ += sqrt(greaterElementFromColumn);
        } else{
            maxZ += log10(greaterElementFromColumn);
        }
        greaterElementFromColumn += chromosomeSize;
    } 

    return maxZ;
}

float fitness(vector<double> individuo, int chromosomeSize, int populationSize){
    
    float rainhasNaoAtacadas, funcaoObjetivo, fit;

    rainhasNaoAtacadas = float(howManyQueensCantBeAttacked(individuo, chromosomeSize));
   
    funcaoObjetivo = funcObjetivo(individuo, chromosomeSize);

    fit = (funcaoObjetivo / maxZToOurContext(chromosomeSize)) + (-1) * fmax(0, (chromosomeSize - rainhasNaoAtacadas) / chromosomeSize);

    if (fit < 0) fit = 0;

    return fit;

}

int find(int valor, vector<double> fellow){
    for (int i = 0; i < int(fellow.size()); i++){
        if(fellow[i] == valor){
            return i; //se tem o valor no vetor, retorna a posição
        }
    }
    return -1; //se nao tem o valor no vetor, retorna -1
}

vector<vector<double>> crossOverPMX(vector<double> fellow1, vector<double> fellow2, int chromosomeSize){

    //double eachChance = (double) 1 / chromosomeSize;
    //double summChance = 0;
    vector<double> newFellow1;
    vector<double> newFellow2;
    vector<double> auxFellow1;
    vector<double> auxFellow2;
    vector<double> auxFellow;
    vector<vector<double>> newFellows;
    int cutIndex, cutIndex2;
    int aux, aux2, valor;
    int cont = 0;

    uniform_real_distribution<double> distribuicaoInteira{0, double(chromosomeSize)};

    //cout << "PMX: Crossover function" << endl;

    //Passo 1: Escolhendo aleatoriamente dois pontos para corte
    cutIndex = distribuicaoInteira(engine);
    cutIndex2 = distribuicaoInteira(engine);

    //if(cutIndex == cutIndex2){
     //   cutIndex2 = distribuicaoInteira(engine);
    //}

    if (cutIndex > cutIndex2){
        double aux = cutIndex;
        cutIndex = cutIndex2;
        cutIndex2 = aux;
    }

    //cout << "cutIndex = " << cutIndex << endl;
    //cout << "cutIndex2 = " << cutIndex2 << endl;

    //Passo 2: Trocar as partes da matching section
    for(int i = 0; i < chromosomeSize; i++){
        if(i < cutIndex){
            newFellow1.push_back(fellow1[i]);
            newFellow2.push_back(fellow2[i]);
        }else if (i == cutIndex || i <= cutIndex2){
            newFellow1.push_back(fellow2[i]);
            newFellow2.push_back(fellow1[i]);

            auxFellow1.push_back(fellow1[i]);
            auxFellow2.push_back(fellow2[i]);
        }
        else{
            newFellow1.push_back(fellow1[i]);
            newFellow2.push_back(fellow2[i]);
        }
    }


    //Passo 3: Mapear o restante dos alelos
    //Verifica se os alelos que estão fora da matching section corresponde 
    //a algum dos alelos dentro da matching section
    //deixei uma imagem no drive pra tentar ajudar a identificar as variaveis
    
    for(int i = 0; i < chromosomeSize; i++){
        if(i < cutIndex || i > cutIndex2){
            
            for(int j = 0; j < 2; j++){ //for para os dois filhos
                
                if(j == 0){ //Para filho 1
                    aux = find(fellow1[i], auxFellow2);
                }else{//Para filho 2
                    aux = find(fellow2[i], auxFellow1);
                }

                if(aux != -1){

                    while(1){
                        if(cont == 0 && j == 0){ //Para filho 1
                            valor = auxFellow1[aux];
                            auxFellow = auxFellow2;
                        }else if(cont != 0 && j == 0){ //Para filho 1
                            valor = auxFellow1[aux2];
                            auxFellow = auxFellow2;
                        }else if(cont == 0 && j == 1){ //Para filho 2
                            valor = auxFellow2[aux];
                            auxFellow = auxFellow1;
                        }else if(cont != 0  && j == 1){ //Para filho 2
                            valor = auxFellow2[aux2];
                            auxFellow = auxFellow1;
                        }

                        aux2 = find(valor, auxFellow);

                        if(aux2 == -1){
                            if(j == 0){ //Para filho 1
                                newFellow1[i] = valor;
                            }else{ //Para filho 2
                                newFellow2[i] = valor;
                            }
                            
                            cont = 0;
                            auxFellow.clear();
                            break;
                        }else{
                            cont += 1;
                            auxFellow.clear();

                        }        
                    }
                }
            }
        }
    }

    newFellows.push_back(newFellow1);
    newFellows.push_back(newFellow2);

    return newFellows;
}

vector<double> mutation(vector<double> fellow, int chromosomeSize){
    vector<double> newFellow = fellow;
    int new_pos;
    int val_atual;
    int val_pos;
    uniform_real_distribution<double> distribuicaoInteira{0, double(chromosomeSize)};

    //cout << "Mutation:" << endl;
    
    for (int i = 0; i < chromosomeSize; i++){
        double aleatorio = distribuicaoBool(engine);
        if (aleatorio <= 0.05){
            
            new_pos = distribuicaoInteira(engine);
            val_pos = newFellow[new_pos];
            val_atual = newFellow[i];

            if(new_pos == i){ //posição deve ser diferente da atual
                new_pos = distribuicaoInteira(engine);
                val_pos = newFellow[new_pos];
            }else{
                //cout << "Posição " << i << " troca com a posição " << new_pos << endl;

                newFellow[new_pos] = val_atual;
                newFellow[i] = val_pos;

            }
            
        }
    }

    return newFellow;
}

vector<int> melhor_pior_individuo(vector<float> fitnessPerIndiv, int populationSize){
    float val_melhor_ind = fitnessPerIndiv[0];
    float val_pior_ind = fitnessPerIndiv[0];
    int posicao_melhor = 0;
    int posicao_pior = 0;
    vector<int> individuos;

    for(int i=0; i<populationSize; i++){
        if(fitnessPerIndiv[i] > val_melhor_ind){
            posicao_melhor = i;
            val_melhor_ind = fitnessPerIndiv[i];
        }

        if(fitnessPerIndiv[i] < val_pior_ind){
            posicao_pior = i;
            val_pior_ind = fitnessPerIndiv[i];
        }
    }

    individuos.push_back(posicao_melhor);
    individuos.push_back(posicao_pior);

    return individuos;
}

vector<vector<double>> selectionRoutineShortExplained(vector<vector<double>> population, int populationSize, int chromosomeSize, vector<float> fitnessPerIndiv_QueensNotAttacked){
    // TO DO: aqui dentro que deve ser chamado o how many queens attacked

    vector<float> probabilities; //probabilidades entre 0 e 1, flutuante.
    float totalSumFitness = 0;
    float localFitnessProportion = 0;
    float actualSumm = 0;
    int indexChoosedToFirstPair = 0;
    int indexChoosedToSecondPair = 0;

    vector<vector<double>> generatePopAfterCrossover;
    vector<vector<double>> auxFellows;

    for (int i = 0; i< populationSize; i++){
        totalSumFitness += fitnessPerIndiv_QueensNotAttacked[i];
    }

    //cout << endl << "Total lucros com todos cromossomos para cálculo da proporção individual do fitness: " << totalSumFitness << endl;
    //cout << "Proporção de fitness por cromossomo (fitness relativo):" << endl;
    //for (int i = 0; i< populationSize; i++){
    //    cout << "Cr[" << i << "]: " << fitnessPerIndiv_QueensNotAttacked[i]/totalSumFitness << endl;
    //}

    // OBS: sempre que i % 2 ==0, então é primeiro do par, logo não precisa verificar e adicionar ao int indexChoosed
    // Se for impar, precisa verificar o indexChoosed e gerar até alcançar o objetivo.

    //cout << "Rotina de seleção dos novos pares de indivíduos/cromossomos" << endl;

    // Neste for se itera pela população para geração dosnovos indivíduos, a partir do crossover
    for (int i = 0; i < populationSize / 2; i++){
       
        //Geração de cada indivíduo a partir de 1 par de indívduos
        for (int j=0; j < 2; j++){
            double aleatorio = distribuicaoBool(engine);
            actualSumm = 0;
            for (int k=0; k<populationSize; k++){ //Itera por cada individuo para avalair sua fitness para sorteio
                localFitnessProportion = fitnessPerIndiv_QueensNotAttacked[k] / totalSumFitness;
                actualSumm += localFitnessProportion;
                if (aleatorio < actualSumm && j == 0){ //Se primeiro elemento do par
                    indexChoosedToFirstPair = k;
                    break;
                }
                else if(aleatorio < actualSumm && j == 1 && k == indexChoosedToFirstPair ){ //Caso segundo elemento do par e repita
                    aleatorio = distribuicaoBool(engine);
                    actualSumm = 0;
                    k = 0;
                }
                else if (aleatorio < actualSumm && j == 1 && k != indexChoosedToFirstPair){ // Caso segundo elemento do par e não é repetição
                    indexChoosedToSecondPair = k;
                    break;
                }
            }
        }
        //cout << endl << "Os indivíduos/cromossomos que dão origem a este par são: os de índice " << indexChoosedToFirstPair << " e " << indexChoosedToSecondPair << endl;

        // AQUI REALIZA PROCESSO, PRIMEIRO CROSSOVER PMX, CASO CAIA NA PROBABILIDADE, DEPOIS MUTAÇÃO EM CADA UM 
        double crossOverChance = distribuicaoBool(engine);
        if (crossOverChance >= 0.2){ // Chance de realizar crossover: 80%
            //cout << "Em it" << i << " vai realizar crossover!" << endl;

            auxFellows = crossOverPMX(population[indexChoosedToFirstPair], population[indexChoosedToSecondPair], chromosomeSize);
            generatePopAfterCrossover.push_back(auxFellows[0]);
            generatePopAfterCrossover.push_back(auxFellows[1]);
            auxFellows.clear();
                
        }else{
            generatePopAfterCrossover.push_back(population[indexChoosedToFirstPair]);
            generatePopAfterCrossover.push_back(population[indexChoosedToSecondPair]);
            //cout << "Em it" << i << " NÃO vai realizar crossover!" << endl;
        }
    
    }

    // IF POPULATION SIZE % 2 == 1, ENTÃO FAZ UM GERAÇÃO EXTRA SÓ COM MUTAÇÃO, SEM CROSSOVER, SOMENTE MUTAÇÃO!
    if(populationSize % 2 == 1){
        //Giro a roleta para elemento único, sem par, enviando-a para mutação direto.
        double aleatorio = distribuicaoBool(engine);
        actualSumm = 0;
        for (int k=0; k<populationSize; k++){ //Itera por cada individuo para avalair sua fitness para sorteio
            localFitnessProportion = fitnessPerIndiv_QueensNotAttacked[k] / totalSumFitness;
            actualSumm += localFitnessProportion;
            if (aleatorio < actualSumm){
                indexChoosedToFirstPair = k;
                break;
            }
        }
        generatePopAfterCrossover.push_back(population[indexChoosedToFirstPair]);
    }

    //printPopulation(generatePopAfterCrossover, populationSize, chromosomeSize, "Printando população pós crossover:");

    //Depois de passar pelo crossover vem a mutação:
    vector<vector<double>> generatePopAfterMutation;
    vector<double> mutFellows;

    for (int i = 0; i < populationSize; i++){
    
        mutFellows = mutation(generatePopAfterCrossover[i], chromosomeSize);
        generatePopAfterMutation.push_back(mutFellows);
        mutFellows.clear();
       
     }

    printPopulation(generatePopAfterMutation, populationSize, chromosomeSize, "População final desta geração:");
    
    return generatePopAfterMutation;


}

vector<vector<double>> selectionRoutineTorneio(vector<vector<double>> population, int populationSize, int chromosomeSize, vector<float> fitnessPerIndiv_QueensNotAttacked){
    float totalSumFitness = 0;
    float localFitnessProportion = 0;
    float actualSumm = 0;
    int indexChoosedToFirstPair = 0;
    int indexChoosedToSecondPair = 0;

    vector<vector<double>> generatePopAfterCrossover;
    vector<vector<double>> auxFellows;

    //Para método do TORNEIO especificamente;
    uniform_real_distribution<double> distribuicaoInteira{0, double(populationSize)};
    vector<int> indexes;
    vector<float> fitness;
    vector<int> m_p_ind;
    int k = 5;
    int aux;
    double kp = 1;


    for (int i = 0; i< populationSize; i++){
        totalSumFitness += fitnessPerIndiv_QueensNotAttacked[i];
    }

    // Neste for se itera pela população para geração dos novos indivíduos, a partir do crossover
    for (int i = 0; i < populationSize / 2; i++){
       
        //Seleção usando Torneio
        for(int j = 0; j < 2; j++){

            //Passo 1: seleciona-se k indivíduos aleatoriamente. Usualmente tem-se k = 2.
            for (int x = 0; x < k; x++){
                aux = distribuicaoInteira(engine);

                indexes.push_back(aux); 
            }

            //Dada a probabilidade kp, escolhe-se o melhor ou pior indivíduo
            for (int x = 0; x < k; x++){  //Pegando a fitness dos k valores
                fitness.push_back(fitnessPerIndiv_QueensNotAttacked[int(indexes[x])]);
            }

            m_p_ind = melhor_pior_individuo(fitness, k); //retorna a posição dentro do vetor fitness**

            if(j == 0){
                if(kp >= distribuicaoBool(engine)){//Escolhe o melhor indivíduo
                    indexChoosedToFirstPair = indexes[m_p_ind[0]];
                }else{//Escolhe o pior indivíduo
                    indexChoosedToFirstPair = indexes[m_p_ind[1]];
                }
            }else{
                if(kp >= distribuicaoBool(engine)){//Escolhe o melhor indivíduo
                    indexChoosedToSecondPair = indexes[m_p_ind[0]];
                }else{//Escolhe o pior indivíduo
                    indexChoosedToSecondPair = indexes[m_p_ind[1]];
                }
            }

            indexes.clear();
            fitness.clear();
            m_p_ind.clear();
            
            
        }

        //cout << "index1 = " << indexChoosedToFirstPair << endl;
        //cout << "index2 = " << indexChoosedToSecondPair << endl;

        // AQUI REALIZA PROCESSO, PRIMEIRO CROSSOVER PMX, CASO CAIA NA PROBABILIDADE, DEPOIS MUTAÇÃO EM CADA UM 
        double crossOverChance = distribuicaoBool(engine);

        if (crossOverChance >= 0.2){ // Chance de realizar crossover: 80%
            //cout << "Em it" << i << " vai realizar crossover!" << endl;

            auxFellows = crossOverPMX(population[indexChoosedToFirstPair], population[indexChoosedToSecondPair], chromosomeSize);
            generatePopAfterCrossover.push_back(auxFellows[0]);
            generatePopAfterCrossover.push_back(auxFellows[1]);
            auxFellows.clear();
                
        }else{
            generatePopAfterCrossover.push_back(population[indexChoosedToFirstPair]);
            generatePopAfterCrossover.push_back(population[indexChoosedToSecondPair]);
            //cout << "Em it" << i << " NÃO vai realizar crossover!" << endl;
        }
    
    }

    // IF POPULATION SIZE % 2 == 1, ENTÃO FAZ UM GERAÇÃO EXTRA SÓ COM MUTAÇÃO, SEM CROSSOVER, SOMENTE MUTAÇÃO!
    if(populationSize % 2 == 1){
        //Giro a roleta para elemento único, sem par, enviando-a para mutação direto.
        double aleatorio = distribuicaoBool(engine);
        actualSumm = 0;
        for (int k=0; k<populationSize; k++){ //Itera por cada individuo para avalair sua fitness para sorteio
            localFitnessProportion = fitnessPerIndiv_QueensNotAttacked[k] / totalSumFitness;
            actualSumm += localFitnessProportion;
            if (aleatorio < actualSumm){
                indexChoosedToFirstPair = k;
                break;
            }
        }
        generatePopAfterCrossover.push_back(population[indexChoosedToFirstPair]);
    }

    //printPopulation(generatePopAfterCrossover, populationSize, chromosomeSize, "Printando população pós crossover:");

    //Depois de passar pelo crossover vem a mutação:
    vector<vector<double>> generatePopAfterMutation;
    vector<double> mutFellows;

    for (int i = 0; i < populationSize; i++){
    
        mutFellows = mutation(generatePopAfterCrossover[i], chromosomeSize);
        generatePopAfterMutation.push_back(mutFellows);
        mutFellows.clear();
       
     }

    printPopulation(generatePopAfterMutation, populationSize, chromosomeSize, "População final desta geração:");
    
    return generatePopAfterMutation;


}

void objectiveFunctionOutput(double sum, vector<double> elements, double qtdElements, string frase){
    //cout << "sum = " << sum << endl;
    //cout << "qtdElements = " << qtdElements << endl;
    double sdPow = 0, sd = 0;
    double mean = sum / qtdElements;

    for (int i = 0; i < qtdElements; i++){
        sdPow += pow(elements[i] - mean, 2);
    }

    sd = sqrt(sdPow / qtdElements);

    cout << frase << mean << " +- " << sd << endl;
    
    cout << "## ## ## ## ## ## ## ##" << endl;
}

void verificaRestricao(int qtdQueensAttacked, int chromosomeSize){
    if(qtdQueensAttacked < chromosomeSize){
        cout << "Essa solução é inválida porque têm rainhas sendo atacadas." << endl;
    }else{
        cout << "Essa solução é válida porque nenhuma rainha está sendo atacada." << endl;
    }
}


int main(int argc, char const *argv[]) {

    // Partindo do princípio que quero uma config de 10 variáveis
    //int quantidade_variaveis = 10;

    if (argc != 6){
		printf("5 argumentos: <COD> <RUN> <GEN> <POP> <DIM> \n");
		return 1;
	}


    int codificationType; // (COD) - tipo de codificação 0 - bin; 1 - int; 2 - int-perm; 3 - real
    int executionsNumber; // (RUN) - número de execuções
    int numberGenerations; // (GEN) - número de gerações
    int populationSize; // (POP) - tamanho da população
    int chromosomeSize; // (DIM) - tamanho do cromossomo

    sscanf(argv[1], "%d", &codificationType);
    sscanf(argv[2], "%d", &executionsNumber);
    sscanf(argv[3], "%d", &numberGenerations);
    sscanf(argv[4], "%d", &populationSize);
    sscanf(argv[5], "%d", &chromosomeSize);

    vector<float> bestFitnessPerGeneration;
    vector<float> meanFitnessPerGeneration;
    vector<vector<float>> bestMeanFitnessPerGeneration;
    vector<vector<float>> meanFitnessPopPerGeneration;
    vector<vector<float>> bestFitnessPerGenerationToAllExecs;
    vector<vector<float>> meanFitnessPerGenerationToAllExecs; 
    vector<float> meanFitnessIndivPerExecution;
    vector<float> meanMeanFitnessPopPerExecution;
    vector<float> bestFitnessPerExec;

    vector<double> bestFellowAllExec;
    float bestFitnessAllExec = 0;

    vector<double> eachObjFunctionValueAllExec;
    double sumObjFunctionValueAllExec = 0;

    vector<double> eachRainhasNaoAtacadasValueAllExec;
    double sumRainhasNaoAtacadasValueAllExec = 0;
    
   // fitness(chromosomeSize, populationSize);

    // Aqui vai um for externo para cada geração possível 
    for (int numExec = 0; numExec < executionsNumber; numExec++){

        vector<vector<double>> pop_inicial = gerarConfiguracaoAleatoria(codificationType, executionsNumber, numberGenerations, populationSize, chromosomeSize);
        vector<vector<double>> pop_atual = pop_inicial;
        vector<vector<double>> pop_atual_aux;

        vector<float> fitnessPerIndiv;

        // Dados do melhor usuário:
        vector<double> bestFellow;
        float bestFitness = 0;

        float sumFitness = 0;
        float sumBestFitnessIndiv = 0;
        float sumAverageFitnessPop = 0;

        double sumObjFunctionValuePerExec = 0;
        vector<double> eachObjFunctionValuePerExec;

        double sumRainhasNaoAtacadasValuePerExec = 0;
        vector<double> eachRainhasNaoAtacadasValuePerExec;

        // População inicial:
        cout << endl;
        cout << "########### EXECUÇÃO " << numExec << ": ########### " << endl;
        printPopulation(pop_inicial, populationSize, chromosomeSize, "Printando população inicial:");

        ofstream generation_by_exec_files;
        string outputt = "RainhasProblem_fitnessOutputFile_" + to_string(chromosomeSize) + "rainhas_" + to_string(numberGenerations) + "gens_" +to_string(executionsNumber) +"execs_ConvExec" + to_string(numExec+1) +".txt";
        generation_by_exec_files.open(outputt);
        generation_by_exec_files << "Geracao MelhorFitn MediaFitness\n";

        for(int generation = 0; generation < numberGenerations; generation ++){
            cout << endl;
            cout << "########### GERAÇÃO " << generation << ": ########### " << endl;

            printPopulation(pop_atual, populationSize, chromosomeSize, "População do início desta geração:");
        
            for (int i = 0; i < populationSize; i++){

                // FITNESS
                float howManyCantBeAttacked = float(howManyQueensCantBeAttacked(pop_atual[i], chromosomeSize));
                float val_fit = fitness(pop_atual[i], chromosomeSize, populationSize);
                float funcObj = funcObjetivo(pop_atual[i], chromosomeSize);
                
                //float val_fit = howManyCantBeAttacked / float(chromosomeSize) ;
                sumFitness += val_fit;
                fitnessPerIndiv.push_back(val_fit);
                cout << "Quantidade de rainhas que não foram atacadas: " << howManyCantBeAttacked << "/" << chromosomeSize << " (" << 
                fitnessPerIndiv[i]<< ")" <<  endl;    

                sumRainhasNaoAtacadasValuePerExec += howManyCantBeAttacked;
                sumRainhasNaoAtacadasValueAllExec += howManyCantBeAttacked;
                eachRainhasNaoAtacadasValuePerExec.push_back(howManyCantBeAttacked);
                eachRainhasNaoAtacadasValueAllExec.push_back(howManyCantBeAttacked);

                sumObjFunctionValuePerExec += funcObj;
                sumObjFunctionValueAllExec += funcObj;
                eachObjFunctionValuePerExec.push_back(funcObj);
                eachObjFunctionValueAllExec.push_back(funcObj);
                
            }

            vector<int> m_p_ind = melhor_pior_individuo(fitnessPerIndiv, populationSize);

            // Caso haja um novo melhor indivíduo histórico, pegamos este novo:
            if (fitnessPerIndiv[m_p_ind[0]] > bestFitness){
                bestFellow.clear();
                for (int i = 0; i < chromosomeSize; i++){  
                    bestFellow.push_back(pop_atual[m_p_ind[0]][i]);
                }
                bestFitness = fitnessPerIndiv[m_p_ind[0]];
                cout << endl;
                printFellow(bestFellow, populationSize, chromosomeSize, "Novo melhor indivíduo: ");
                cout << "Sua fitness: " << bestFitness << endl;

                // Se o melhor dessa geração for um novo melhor gerado
                bestFitnessPerGeneration.push_back(fitnessPerIndiv[m_p_ind[0]]);
                sumBestFitnessIndiv += fitnessPerIndiv[m_p_ind[0]];

                // Se o novo melhor indivíduo é melhor que o melhor de todas execs:
                if (bestFitness > bestFitnessAllExec){
                    bestFellowAllExec.clear();
                    for (int i = 0; i < chromosomeSize; i++){  
                        bestFellowAllExec.push_back(bestFellow[i]);
                    }
                    bestFitnessAllExec = bestFitness;
                }
            }
            else{
                // Caso nessa geração não haja um novo melhor, o melhor histórico substitui o pior, para garantir sua manutenção

                for (int i = 0; i < chromosomeSize; i++){ // corrigir o cromossomo do pior indivíduo
                    pop_atual[m_p_ind[1]][i] = bestFellow[i];
                }
                sumFitness = sumFitness - fitnessPerIndiv[m_p_ind[1]] + bestFitness; // corrigir sumFitness
                fitnessPerIndiv[m_p_ind[1]] = bestFitness; // corrigir fitness do pior indivíduo

                // Se o melhor dessa geração for o pior substituído pelo melhor:
                bestFitnessPerGeneration.push_back(fitnessPerIndiv[m_p_ind[1]]);
                sumBestFitnessIndiv += fitnessPerIndiv[m_p_ind[1]];
            }

            generation_by_exec_files << to_string(generation + 1) + " " + to_string(bestFitness) + " " + to_string(sumFitness / populationSize) + "\n";


            //bestFitnessPerGeneration.push_back(fitnessPerIndiv[m_p_ind[0]]);
            //sumBestFitnessIndiv += fitnessPerIndiv[m_p_ind[0]];
            meanFitnessPerGeneration.push_back(sumFitness / populationSize);
            sumAverageFitnessPop += sumFitness / populationSize;
            sumFitness = 0;

            cout << endl;
            cout << "O melhor indivíduo é o Cr.["<< m_p_ind[0] <<"] com fitness = " << fitnessPerIndiv[m_p_ind[0]]<< endl;
            
            //Esse selection vai precisar retorna vector<vector<double>> e este será recebido no final de cada geração para usar como geração base seguinte para próximo loop.
            pop_atual_aux = selectionRoutineTorneio(pop_atual, populationSize, chromosomeSize, fitnessPerIndiv);
            //pop_atual_aux = selectionRoutineShortExplained(pop_atual, populationSize, chromosomeSize, fitnessPerIndiv);
            
            pop_atual.clear();
            pop_atual = pop_atual_aux;

            pop_atual_aux.clear();
            fitnessPerIndiv.clear();
            m_p_ind.clear();
        }

        generation_by_exec_files.close();

        // Insere dados no objetos que vai gerar o gráfico de convergência médio de todas execuções (novo):
        //bestMeanFitnessPerGeneration.push_back(bestFitnessPerGeneration);
        //meanFitnessPopPerGeneration.push_back(meanFitnessPopPerGenerationSum);

        meanFitnessIndivPerExecution.push_back(sumBestFitnessIndiv / numberGenerations);
        meanMeanFitnessPopPerExecution.push_back(sumAverageFitnessPop / numberGenerations);

        bestFitnessPerGenerationToAllExecs.push_back(bestFitnessPerGeneration);
        meanFitnessPerGenerationToAllExecs.push_back(meanFitnessPerGeneration);

        bestFitnessPerGeneration.clear();
        meanFitnessPerGeneration.clear();

        // EXIBIÇÃO MELHOR INDIVÍDUO ENCONTRADO NESTA EXEC:
        printFellow(bestFellow, populationSize, chromosomeSize, "### MELHOR INDIVÍDUO ENCONTRADO NESTA EXECUÇÃO: ### ");
        //cout << fitnessVerbose(bestFellow, chromosomeSize) << "/1 de fitness" << endl;
        //int fobj = howManyQueensCantBeAttacked(bestFellow, chromosomeSize);
        float fit = fitness(bestFellow, chromosomeSize, populationSize);
        float fobj = funcObjetivo(bestFellow, chromosomeSize);
        int howManyCantBeAttacked = howManyQueensCantBeAttacked(bestFellow, chromosomeSize);
        //double fit = double(fobj) / double(chromosomeSize);
        cout << "Fitness: " << fit << endl;
        cout << "Função objetivo (Z): " << fobj << "/" << maxZToOurContext(chromosomeSize) << endl;
        cout << "Quantidade rainhas não atacadas: " << howManyCantBeAttacked << "/" << chromosomeSize << endl;
        verificaRestricao(howManyCantBeAttacked, chromosomeSize);

        bestFitnessPerExec.push_back(fit);

        // EXIBIÇÃO DA MÉDIA E DESVIO PADRÃO DA FUNÇÃO OBJETIVO:
        cout << endl << "## ESTATÍSTICAS DA FUNÇÃO OBJETIVO DESTA EXECUÇÃO: ##" << endl;
        objectiveFunctionOutput(sumObjFunctionValuePerExec, eachObjFunctionValuePerExec, populationSize * numberGenerations, "Média da função objetivo com desvio padrão: ");
        
        // EXIBIÇÃO DA MÉDIA E DESVIO PADRÃO DA QTDD DE RAINHAS NÃO ATACADAS:
        cout << endl << "## ESTATÍSTICAS DAS RAINHAS NÃO ATACADAS DESTA EXECUÇÃO: ##" << endl;
        objectiveFunctionOutput(sumRainhasNaoAtacadasValuePerExec, eachRainhasNaoAtacadasValuePerExec, populationSize * numberGenerations, "Média da qtdd de rainhas não atacadas com desvio padrão: ");
        
        
    }

    // EXIBIÇÃO MELHOR INDIVÍDUO ENCONTRADO TODAS EXECS:
    cout << endl;
    cout << "################### GERAL ###################" << endl;
    cout << endl;
    printFellow(bestFellowAllExec, populationSize, chromosomeSize, "### MELHOR INDIVÍDUO ENCONTRADO EM TODAS AS EXECUÇÕES: ### ");
    //cout << fitnessVerbose(bestFellowAllExec, chromosomeSize) << "/1 de fitness " << endl;
    float fitAllExec = fitness(bestFellowAllExec, chromosomeSize, populationSize);
    float fobjAllExec = funcObjetivo(bestFellowAllExec, chromosomeSize);
    int howManyCantBeAttackedAllExec = howManyQueensCantBeAttacked(bestFellowAllExec, chromosomeSize);
    cout << "Fitness: " << fitAllExec << endl;
    cout << "Função objetivo (Z): " << fobjAllExec << "/" << maxZToOurContext(chromosomeSize) << endl;
    cout << "Quantidade rainhas não atacadas: " << howManyCantBeAttackedAllExec << "/" << chromosomeSize << endl;
    verificaRestricao(howManyCantBeAttackedAllExec, chromosomeSize);

    // EXIBIÇÃO DA MÉDIA E DESVIO PADRÃO DA FUNÇÃO OBJETIVO:
    cout << endl << "## ESTATÍSTICAS DA FUNÇÃO OBJETIVO PARA TODAS EXECUÇÕES: ##" << endl;
    objectiveFunctionOutput(sumObjFunctionValueAllExec, eachObjFunctionValueAllExec, populationSize * numberGenerations * executionsNumber, "Média da função objetivo com desvio padrão: ");

    // EXIBIÇÃO DA MÉDIA E DESVIO PADRÃO DA QTDD DE RAINHAS NÃO ATACADAS:
    cout << endl << "## ESTATÍSTICAS DA QTDD DE RAINHAS NÃO ATACADAS PARA TODAS EXECUÇÕES: ##" << endl;
    objectiveFunctionOutput(sumRainhasNaoAtacadasValueAllExec, eachRainhasNaoAtacadasValueAllExec, populationSize * numberGenerations * executionsNumber, "Média da qtdd de rainhas não atacadas com desvio padrão: ");
   

    // OBS: caso seja apenas uma execução, tirei pois lá em cima cada execução eu já gero arq mesmo por geração

    // Caso mais de uma  execução, exibe e grava em arquivo a média do fitness de cada exec e a média da média do fitness da pop de cada execução :
    if(executionsNumber > 1){
        ofstream fitness_file;
        string output = "RainhasProblem_fitnessOutputFile_" + to_string(chromosomeSize) + "rainhas_" + to_string(numberGenerations) + "gens_" +to_string(executionsNumber) +"execs.txt";
        fitness_file.open(output);

        cout << endl << "LISTA DE CONVERGÊNCIA POR EXECUÇÃO" << endl;
        cout << "Execução | Média melhor fitn.  | Média média fitness pop | Melhor Fit Exec" << endl;
        fitness_file << "Execucao MediaMelhorFitn MediaMediaFitnessPop MelhorFit\n";
        

        for (int exec = 0; exec < executionsNumber; exec++){
            cout << exec+1 << "\t | " << meanFitnessIndivPerExecution[exec] << "\t | " << meanMeanFitnessPopPerExecution[exec] << "\t | "  << bestFitnessPerExec[exec] << endl; 
            fitness_file << to_string(exec+1) + " " + to_string(meanFitnessIndivPerExecution[exec]) + " " + to_string(meanMeanFitnessPopPerExecution[exec]) + " " + to_string(bestFitnessPerExec[exec]) + "\n";
            
        }

        fitness_file.close();

    }

    // GERAÇÃO DE GRÁFICO DE CONVERGÊNCIA DE FITNESS COM MÉDIA DA GERAÇÃO DE CADA EXECUÇÃO (PARA 1 OU MAIS EXECUÇÕES)
    ofstream convergenceToAllExec_file;
    string output_convergence = "RainhasProblem_fitnessOutputFile_" + to_string(chromosomeSize) + "rainhas_" + to_string(numberGenerations) + "gens_" +to_string(executionsNumber) +"execs_ConvergenceMeanAllExecs.txt";
    convergenceToAllExec_file.open(output_convergence);
    convergenceToAllExec_file << "Geracao MelhorFitn MediaFitness\n";

    for (int ger = 0; ger < numberGenerations; ger++){

        float sumBestFitness = 0;
        float sumMeanFitnessPop = 0;
        float meanBestFitnessAllExecs = 0;
        float meanMeanFitnessPopAllExecs = 0;

        for (int exec = 0; exec < executionsNumber; exec++){
            sumBestFitness += bestFitnessPerGenerationToAllExecs[exec][ger];
            sumMeanFitnessPop += meanFitnessPerGenerationToAllExecs[exec][ger];
        }

        meanBestFitnessAllExecs = sumBestFitness / 10;
        meanMeanFitnessPopAllExecs = sumMeanFitnessPop / 10;

        convergenceToAllExec_file << to_string(ger + 1) + " " + to_string(meanBestFitnessAllExecs) + " " + to_string(meanMeanFitnessPopAllExecs) + "\n";

        sumBestFitness = 0;
        sumMeanFitnessPop = 0;
        
    }

    convergenceToAllExec_file.close();

    

}
