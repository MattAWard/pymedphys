// cypress.config.js
//
// Migration from cypress.json (Cypress 9) to cypress.config.js (Cypress 15)
//
// Key changes:
// 1. cypress.json is no longer supported from Cypress 10+
// 2. Test directory moved: cypress/integration/ -> cypress/e2e/
// 3. Spec file pattern: **/*.cy.js (rename from *.js is recommended but optional)
// 4. Plugins file (cypress/plugins/index.js) replaced by setupNodeEvents below
// 5. cypress-file-upload plugin replaced by built-in cy.selectFile()
//
const { defineConfig } = require("cypress");

module.exports = defineConfig({
  // Carried over from old cypress.json
  viewportWidth: 1000,
  viewportHeight: 1200,
  projectId: "tgt8f6",

  // Increase default command timeout for Streamlit app loading
  defaultCommandTimeout: 30000,

  // Recommended: disable video for faster CI runs (enable if needed for debugging)
  video: false,

  // Retry configuration for flaky Streamlit tests
  retries: {
    runMode: 2,
    openMode: 0,
  },

  e2e: {
    // Base URL for the Streamlit app
    baseUrl: "http://localhost:8501",

    // Test file pattern — files moved from cypress/integration/ to cypress/e2e/
    // Using .cy.js extension per Cypress 10+ convention
    specPattern: "cypress/e2e/**/*.cy.js",

    // Support file — this replaces the old cypress/support/index.js
    supportFile: "cypress/support/e2e.js",

    setupNodeEvents(on, config) {
      // This replaces cypress/plugins/index.js
      // Add any event listeners or task definitions here
      // Example:
      // on('task', { log(message) { console.log(message); return null; } })
      return config;
    },
  },
});
