let geneData = [];

// Load gene data when the page loads
window.onload = async function() {
    try {
        const response = await fetch('gene_data.json');
        geneData = await response.json();
    } catch (err) {
        console.error("Error fetching gene data:", err);
    }
}
function addRow() {
    const tbody = document.getElementById("inputTable").querySelector("tbody");
    const row = tbody.insertRow();
    const cell1 = row.insertCell(0);
    const cell2 = row.insertCell(1);
    const cell3 = row.insertCell(2);
    const cell4 = row.insertCell(3);

    cell1.innerHTML = '<input type="text" class="gene-name" placeholder="Gene...">';
    cell2.innerHTML = '<input type="number" class="gene-length" placeholder="1.5" readonly>';
    cell3.innerHTML = '<input type="number" class="wt" />';
    cell4.innerHTML = '<input type="number" class="mutant" />';
}

function loadExample() {
    const geneNames = document.querySelectorAll(".gene-name");
    const genelen = document.querySelectorAll(".gene-length");
    const wts = document.querySelectorAll(".wt");
    const mutants = document.querySelectorAll(".mutant");
    for(let i = 0; i < geneNames.length; i++) {
    wts[i].value = 500;
    mutants[i].value = 500;
    }
}

function calculateTPM() {
    const geneNames = [...document.querySelectorAll(".gene-name")].map(input => input.value);
    const geneLengths = [...document.querySelectorAll(".gene-length")].map(input => Number(input.value));
    const wtCounts = [...document.querySelectorAll(".wt")].map(input => Number(input.value));
    const mutantCounts = [...document.querySelectorAll(".mutant")].map(input => Number(input.value));

    let wtTPMs = calculateTPMForCondition(wtCounts, geneLengths);
    let mutantTPMs = calculateTPMForCondition(mutantCounts, geneLengths);

    displayResults(geneNames, wtTPMs, mutantTPMs);
}

function calculateTPMForCondition(counts, lengths) {
    const rpk = counts.map((count, i) => count / lengths[i]);
    const scalingFactor = rpk.reduce((acc, curr) => acc + curr, 0);
    return rpk.map(value => (value / scalingFactor) * 1000000);
}

function displayResults(geneNames, wtResults, mutantResults) {
    const tbody = document.getElementById("resultsTable").querySelector("tbody");
    tbody.innerHTML = "";  // Clear previous results

    for(let i = 0; i < geneNames.length; i++) {
        const row = tbody.insertRow();
        row.insertCell(0).textContent = geneNames[i];
        row.insertCell(1).textContent = wtResults[i].toFixed(2);
        row.insertCell(2).textContent = mutantResults[i].toFixed(2);
    }
}

// Add listener to gene name inputs
document.getElementById("inputTable").addEventListener('input', function(event) {
    if (event.target && event.target.classList.contains('gene-name')) {
        updateGeneLength(event.target);
    }
});

function updateGeneLength(inputElement) {
    const geneName = inputElement.value;
    const matchedGene = geneData.find(g => g.name === geneName);
    if (matchedGene) {
        const lengthInput = inputElement.closest('tr').querySelector('.gene-length');
        lengthInput.value = matchedGene.length;
    }
}