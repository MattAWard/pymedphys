// cypress/support/e2e.js
//
// Replaces: cypress/support/index.js
//
// In Cypress 10+, the support file for E2E tests must be named e2e.js (or e2e.ts)
// and placed in cypress/support/. The old index.js is no longer auto-loaded.
//

// Import custom commands
import './commands'

// Preserve the original uncaught:exception handler from the old index.js
Cypress.on('uncaught:exception', (err, runnable) => {
  // returning false here prevents Cypress from failing the test
  return false
})
