import { GoogleGenAI } from "@google/genai";
import { ProteinData } from "../types";
import { SPECIES_CONFIG } from "../constants";

const apiKey = process.env.API_KEY || '';
const ai = new GoogleGenAI({ apiKey });

const MODEL_NAME = "gemini-2.5-flash";

export const analyzeProteins = async (
  proteinData: ProteinData[],
  context: string = ""
): Promise<string> => {
  if (!apiKey) {
    return "Error: API Key is missing. Please check your environment configuration.";
  }

  // Calculate quick stats for the AI
  const speciesStats = Object.keys(SPECIES_CONFIG).map(s => {
      const subset = proteinData.filter(p => p.species === s);
      const median = subset.length ? subset.reduce((acc, curr) => acc + curr.log2FoldChange, 0) / subset.length : 0;
      return `${s}: Observed Median Log2FC = ${median.toFixed(2)} (Expected ${SPECIES_CONFIG[s as keyof typeof SPECIES_CONFIG].expectedLog2FC})`;
  }).join("\n");

  const prompt = `
    You are an expert computational biologist and proteomics quality control specialist using LFQ-Bench.
    
    Task: Analyze the following mixed proteome dataset (Human, Yeast, E. coli).
    The goal is to verify if the mixing ratios are accurate compared to the theoretical ground truth.
    
    Experimental Context: "${context}"
    
    Data Summary:
    ${speciesStats}
    
    Please provide a technical quality control report covering:
    1. **Ratio Accuracy**: Are the observed median Log2 Fold Changes close to the expected values for each species? 
       - Human should be ~${SPECIES_CONFIG['Human'].expectedLog2FC}
       - Yeast should be ~${SPECIES_CONFIG['Yeast'].expectedLog2FC}
       - E.coli should be ~${SPECIES_CONFIG['E.coli'].expectedLog2FC}
    2. **Quantification Performance**: Comment on any potential compression of ratios or systematic shifts.
    3. **Recommendations**: If values are off, suggest potential causes (e.g., mixing error, normalization issues, interference).
    
    Format the output in clean Markdown.
  `;

  try {
    const response = await ai.models.generateContent({
      model: MODEL_NAME,
      contents: prompt,
      config: {
        temperature: 0.3, 
      }
    });
    return response.text || "No analysis generated.";
  } catch (error) {
    console.error("Gemini Analysis Error:", error);
    return "Failed to generate analysis. Please try again later or check your API key.";
  }
};

export const chatWithData = async (
  history: { role: 'user' | 'model', text: string }[],
  newMessage: string,
  dataContext: string
): Promise<string> => {
    if (!apiKey) return "API Key missing.";

    try {
        const systemInstruction = `You are a helper for LFQ-Bench, a proteomics benchmarking tool. 
        You have access to a dataset containing a mix of Human, Yeast, and E. coli proteins at specific ratios.
        Answer questions about the data quality, mixing ratios, and general proteomics normalization concepts.
        Current Data Stats: ${dataContext}`;

        const contents = [
            { role: 'user', parts: [{ text: `System Context: ${systemInstruction}` }] },
            ...history.map(h => ({ role: h.role, parts: [{ text: h.text }] })),
            { role: 'user', parts: [{ text: newMessage }] }
        ];

        const response = await ai.models.generateContent({
            model: MODEL_NAME,
            contents: contents as any,
        });

        return response.text || "I couldn't process that.";
    } catch (error) {
        console.error("Chat Error:", error);
        return "Sorry, I encountered an error responding to that.";
    }
}
